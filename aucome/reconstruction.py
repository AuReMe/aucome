#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
usage:
    aucome reconstruction --run=ID [--cpu=INT] [-v]

options:
    --run=ID    Pathname to the comparison workspace.
    --cpu=INT     Number of cpu to use for the multiprocessing (if none use 1 cpu).
    -v     Verbose.
"""

import configparser
import csv
import docopt
import eventlet
import mpwt
import os
import re
import subprocess
import time

from aucome.utils import parse_config_file
from multiprocessing import Pool


def command_help():
    print(docopt.docopt(__doc__))


def reconstruction_parse_args(command_args):
    args = docopt.docopt(__doc__, argv=command_args)
    run_id = args['--run']
    verbose = args['-v']

    if args["--cpu"]:
        nb_cpu_to_use = int(args["--cpu"])
    else:
        nb_cpu_to_use = 1

    run_reconstruction(run_id, nb_cpu_to_use, verbose)


def run_reconstruction(run_id, nb_cpu_to_use, verbose):
    config_data = parse_config_file(run_id)

    pgdb_from_annotation_path = config_data['pgdb_from_annotation_path']
    studied_organisms_path = config_data['studied_organisms_path']
    log_path = config_data['log_path']
    #check for each study if exist PGDB folder in PGDBs folder, if missing RUN ptools
    chronoDepart = time.time()

    mpwt.multiprocess_pwt(input_folder=config_data['studied_organisms_path'],
                            output_folder=pgdb_from_annotation_path,
                            patho_inference=True,
                            dat_creation=True,
                            dat_extraction=True,
                            number_cpu=nb_cpu_to_use,
                            patho_log=log_path,
                            verbose=verbose)

    chrono = (time.time() - chronoDepart)
    partie_entiere, partie_decimale = str(chrono).split('.')
    chrono = ".".join([partie_entiere, partie_decimale[:3]])

    if os.listdir(pgdb_from_annotation_path) == []:
        print('Pathway-Tools inference failed!')
        return
    if verbose:
        print("Pathway-Tools done in: %ss" %chrono)

    create_padmet_sbml_from_pgdb(run_id, nb_cpu_to_use, verbose)


def create_padmet_sbml_from_pgdb(run_id, nb_cpu_to_use, verbose):
    config_data = parse_config_file(run_id)

    padmet_from_annotation_path = config_data['padmet_from_annotation_path']
    study_from_annot_prefix = config_data['study_from_annot_prefix']
    sbml_from_annotation_path = config_data['sbml_from_annotation_path']
    padmet_utils_path = config_data['padmet_utils_path']
    database_path = config_data['database_path']
    pgdb_from_annotation_path = config_data['pgdb_from_annotation_path']

    aucome_pool = Pool(nb_cpu_to_use)

    all_study_name = set(next(os.walk(config_data['studied_organisms_path']))[1])

    all_study_pgdb = dict([(study_name, "{0}/{1}".format(pgdb_from_annotation_path, study_name))
                        if os.path.isdir("{0}/{1}".format(pgdb_from_annotation_path, study_name))
                        else (study_name, '')
                        for study_name in all_study_name])

    study_padmet_data = []
    for study_name in all_study_name:
        padmet_file = "{0}/{1}{2}.padmet".format(padmet_from_annotation_path, study_from_annot_prefix, study_name)
        pgdb_folder = all_study_pgdb[study_name]
        tmp_padmet_data = {'study_name': study_name, 'pgdb_folder': pgdb_folder, 'padmet_utils_path': padmet_utils_path,
                            'verbose': verbose, 'padmet_file': padmet_file, 'database_path': database_path}
        study_padmet_data.append(tmp_padmet_data)
    aucome_pool.map(create_padmet_from_pgdb, study_padmet_data)

    all_study_padmet = dict([(study_name, "{0}/{1}{2}.padmet".format(padmet_from_annotation_path, study_from_annot_prefix, study_name))
                            if os.path.isfile("{0}/{1}{2}.padmet".format(padmet_from_annotation_path, study_from_annot_prefix, study_name))
                            else (study_name, '')
                            for study_name in all_study_name])

    study_sbml_data = []
    for study_name in all_study_padmet:
        sbml_file = "{0}/{1}{2}.sbml".format(sbml_from_annotation_path, study_from_annot_prefix, study_name)
        padmet_file = all_study_padmet[study_name]
        tmp_sbml_data = {'sbml_file': sbml_file, 'padmet_file': padmet_file, 'padmet_utils_path': padmet_utils_path,
                         'study_name': study_name, 'verbose': verbose}
        study_sbml_data.append(tmp_sbml_data)
    aucome_pool.map(create_sbml, study_sbml_data)


def create_padmet_from_pgdb(tmp_padmet_data):
    study_name = tmp_padmet_data['study_name']
    pgdb_folder = tmp_padmet_data['pgdb_folder']
    verbose = tmp_padmet_data['verbose']
    padmet_file = tmp_padmet_data['padmet_file']
    padmet_utils_path = tmp_padmet_data['padmet_utils_path']
    database_path = tmp_padmet_data['database_path']

    if not os.path.isfile(padmet_file) and pgdb_folder:
        if verbose:
            print("Creating padmet from pgdb for %s" %study_name)
        cmds = ["python3",  padmet_utils_path + "/padmet_utils/connection/pgdb_to_padmet.py", "--pgdb", pgdb_folder, "--output", padmet_file,
                "--padmetRef", database_path, "--source=genome", "--extract-gene", "--no-orphan"]

        if verbose:
            cmds.append('-v')

        subprocess.call(cmds)


def create_sbml(tmp_sbml_data):
    sbml_file = tmp_sbml_data['sbml_file']
    padmet_file = tmp_sbml_data['padmet_file']
    study_name = tmp_sbml_data['study_name']
    verbose = tmp_sbml_data['verbose']
    padmet_utils_path = tmp_sbml_data['padmet_utils_path']

    if not os.path.isfile(sbml_file) and padmet_file:
        if verbose:
            print("Creating sbml from padmet for %s" %study_name)
        cmds = ["python3", padmet_utils_path + "/padmet_utils/connection/sbmlGenerator.py", "--padmet", padmet_file,
                "--output", sbml_file, "--sbml_lvl", "3"]

        if verbose:
            cmds.append('-v')

        subprocess.call(cmds)
