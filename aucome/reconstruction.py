#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
usage:
    aucome reconstruction --run=ID [--cpu=INT] [-v] [--vv]

options:
    --run=ID    Pathname to the comparison workspace.
    --cpu=INT     Number of cpu to use for the multiprocessing (if none use 1 cpu).
    -v     Verbose.
    --vv    Very verbose.
"""

import docopt
import logging
import mpwt
import os
import time

from padmet.utils.connection import pgdb_to_padmet, sbmlGenerator

from aucome.utils import parse_config_file
from multiprocessing import Pool

logger = logging.getLogger('aucome')
logger.setLevel(logging.CRITICAL)
logging.getLogger("mpwt").setLevel(logging.CRITICAL)


def command_help():
    print(docopt.docopt(__doc__))


def reconstruction_parse_args(command_args):
    args = docopt.docopt(__doc__, argv=command_args)
    run_id = args['--run']
    verbose = args['-v']
    veryverbose = args['--vv']

    if args["--cpu"]:
        nb_cpu_to_use = int(args["--cpu"])
    else:
        nb_cpu_to_use = 1

    if veryverbose and not verbose:
        verbose = veryverbose

    run_reconstruction(run_id, nb_cpu_to_use, verbose, veryverbose)


def run_reconstruction(run_id, nb_cpu_to_use, verbose, veryverbose=None):
    if verbose:
        logger.setLevel(logging.DEBUG)
        logging.getLogger("mpwt").setLevel(logging.DEBUG)
        print('--- Running reconstruction step ---')
    start_time = time.time()
    config_data = parse_config_file(run_id)

    pgdb_from_annotation_path = config_data['pgdb_from_annotation_path']
    studied_organisms_path = config_data['studied_organisms_path']
    log_path = config_data['log_path']


    taxon_file = None
    if 'taxon_id.tsv' in set(next(os.walk(config_data['studied_organisms_path']))[2]):
        taxon_file = True

    mpwt.multiprocess_pwt(input_folder=studied_organisms_path,
                            output_folder=pgdb_from_annotation_path,
                            patho_inference=True,
                            flat_creation=True,
                            dat_extraction=True,
                            number_cpu=nb_cpu_to_use,
                            patho_log=log_path,
                            taxon_file=taxon_file,
                            verbose=verbose)

    if os.listdir(pgdb_from_annotation_path) == []:
        print('Pathway-Tools inference failed!')
        return

    create_padmet_sbml_from_pgdb(run_id, nb_cpu_to_use, verbose, veryverbose)

    end_time = (time.time() - start_time)
    integer_part, decimal_part = str(end_time).split('.')
    reconstruction_time = ".".join([integer_part, decimal_part[:3]])

    if verbose:
        print("--- reconstruction step done in: %ss ---" %reconstruction_time)

def create_padmet_sbml_from_pgdb(run_id, nb_cpu_to_use, verbose, veryverbose):
    config_data = parse_config_file(run_id)

    padmet_from_annotation_path = config_data['padmet_from_annotation_path']
    study_from_annot_prefix = config_data['study_from_annot_prefix']
    sbml_from_annotation_path = config_data['sbml_from_annotation_path']
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
        tmp_padmet_data = {'study_name': study_name, 'pgdb_folder': pgdb_folder, 'veryverbose': veryverbose,
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
        tmp_sbml_data = {'sbml_file': sbml_file, 'padmet_file': padmet_file,
                         'study_name': study_name, 'verbose': verbose, 'veryverbose': veryverbose}
        study_sbml_data.append(tmp_sbml_data)
    aucome_pool.map(create_sbml, study_sbml_data)


def create_padmet_from_pgdb(tmp_padmet_data):
    study_name = tmp_padmet_data['study_name']
    pgdb_folder = tmp_padmet_data['pgdb_folder']
    verbose = tmp_padmet_data['verbose']
    padmet_file = tmp_padmet_data['padmet_file']
    database_path = tmp_padmet_data['database_path']
    veryverbose = tmp_padmet_data['veryverbose']

    if not os.path.isfile(padmet_file) and pgdb_folder:
        if verbose:
            print("Creating padmet from pgdb for %s" %study_name)
        pgdb_to_padmet.from_pgdb_to_padmet(pgdb_folder=pgdb_folder, padmetRef_file=database_path, source="genome", extract_gene=True, no_orphan=True, verbose=veryverbose, output_file=padmet_file)


def create_sbml(tmp_sbml_data):
    sbml_file = tmp_sbml_data['sbml_file']
    padmet_file = tmp_sbml_data['padmet_file']
    study_name = tmp_sbml_data['study_name']
    verbose = tmp_sbml_data['verbose']
    veryverbose = tmp_sbml_data['veryverbose']

    if not os.path.isfile(sbml_file) and padmet_file:
        if verbose:
            print("Creating sbml from padmet for %s" %study_name)
        sbmlGenerator.padmet_to_sbml(padmet=padmet_file, output=sbml_file, sbml_lvl=3, verbose=veryverbose)
