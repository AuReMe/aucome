#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
usage:
    aucome draft --run=ID [--cpu=INT] [-v]

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


def draft_parse_args(command_args):
    args = docopt.docopt(__doc__, argv=command_args)
    run_id = args['--run']
    verbose = args['-v']

    if args["--cpu"]:
        nb_cpu_to_use = int(args["--cpu"])
    else:
        nb_cpu_to_use = 1

    run_draft(run_id, nb_cpu_to_use, verbose)


def run_draft(run_id, nb_cpu_to_use, verbose):

    aucome_pool = Pool(nb_cpu_to_use)

    config_data = parse_config_file(run_id)

    studied_organisms_path = config_data['studied_organisms_path']
    padmet_from_annotation_path = config_data['padmet_from_annotation_path']
    study_from_annot_prefix = config_data['study_from_annot_prefix']
    networks_path = config_data['networks_path']
    orthology_based_path = config_data['orthology_based_path']
    padmet_utils_path = config_data['padmet_utils_path']
    database_path = config_data['database_path']

    all_study_name = set(next(os.walk(studied_organisms_path))[1])

    all_study_padmet = dict([(study_name, "{0}/{1}{2}.padmet".format(padmet_from_annotation_path, study_from_annot_prefix, study_name))
                          if os.path.isfile("{0}/{1}{2}.padmet".format(padmet_from_annotation_path, study_from_annot_prefix, study_name))
                          else (study_name, '')
                          for study_name in all_study_name])

    study_draft_data = []
    for study_name in all_study_name:
        tmp_study_data = {'study_name': study_name, 'study_padmet': all_study_padmet[study_name], 'networks_path': networks_path,
                            'orthology_based_path': orthology_based_path, 'padmet_utils_path': padmet_utils_path, 'database_path': database_path,
                            'verbose': verbose}
        study_draft_data.append(tmp_study_data)
    aucome_pool.map(create_draft, study_draft_data)


def create_draft(tmp_study_data):
    study_name = tmp_study_data['study_name']
    study_padmet = tmp_study_data['study_padmet']
    verbose = tmp_study_data['verbose']
    networks_path = tmp_study_data['networks_path']
    orthology_based_path = tmp_study_data['orthology_based_path']
    padmet_utils_path = tmp_study_data['padmet_utils_path']
    database_path = tmp_study_data['database_path']

    output = "{0}/{1}.padmet".format(networks_path, study_name)
    if os.path.exists(output):
        if verbose:
            print("%s already exist, skip" %os.path.basename(output))
            return
    else:
        ortho_sbml_folder = "{0}/{1}".format(orthology_based_path, study_name)
        source_tool = "ORTHOFINDER"
        source_category = "ORTHOLOGY"
        if verbose:
            print("Creating %s" %os.path.basename(output))
        if os.path.exists(study_padmet):
            if verbose:
                print("\tStarting from %s" %os.path.basename(study_padmet))
            padmet_path = study_padmet
            if os.path.exists(ortho_sbml_folder):
                cmds = ["python3",  padmet_utils_path + "/padmet_utils/connection/sbml_to_padmet.py", "--padmetRef", database_path, "--sbml", ortho_sbml_folder,
                        "--padmetSpec", padmet_path, "--output", output, "--source_tool", source_tool, "--source_category", source_category]

                if verbose:
                    cmds.append('-v')
            else:
                if verbose:
                    print("\tNo orthology folder.")
                    print(("\tMove {0} in {1}".format(study_name, output)))
                subprocess.call(["cp", padmet_path, output])
                return
        else:
            if verbose:
                print("\tStarting from an empty PADMET")
            cmds = ["python3",  padmet_utils_path + "/padmet_utils/connection/sbml_to_padmet.py", "--padmetRef", database_path, "--sbml", ortho_sbml_folder,
                    "--padmetSpec", output, "--source_tool", source_tool, "--source_category", source_category]
            if verbose:
                cmds.append('-v')
        if os.path.exists(ortho_sbml_folder) and next(os.walk(ortho_sbml_folder))[2]:
            subprocess.call(cmds)
        else:
            if verbose:
                print("\t%s's folder is empty" %study_name)
            return
