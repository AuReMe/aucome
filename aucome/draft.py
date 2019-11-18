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

import docopt
import os
from shutil import copyfile
from padmet.utils.connection import sbml_to_padmet, sbmlGenerator

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
    database_path = config_data['database_path']
    padmet_from_networks_path = config_data['padmet_from_networks_path']
    sbml_from_networks_path = config_data['sbml_from_networks_path']

    all_study_name = set(next(os.walk(studied_organisms_path))[1])

    all_study_padmet = dict([(study_name, "{0}/{1}{2}.padmet".format(padmet_from_annotation_path, study_from_annot_prefix, study_name))
                          if os.path.isfile("{0}/{1}{2}.padmet".format(padmet_from_annotation_path, study_from_annot_prefix, study_name))
                          else (study_name, '')
                          for study_name in all_study_name])

    study_draft_data = []
    for study_name in all_study_name:
        tmp_study_data = {'study_name': study_name, 'study_padmet': all_study_padmet[study_name], 'networks_path': networks_path,
                            'orthology_based_path': orthology_based_path, 'database_path': database_path,
                            'padmet_from_networks_path': padmet_from_networks_path, 'sbml_from_networks_path': sbml_from_networks_path, 'verbose': verbose}
        study_draft_data.append(tmp_study_data)
    aucome_pool.map(create_draft, study_draft_data)


def create_draft(tmp_study_data):
    study_name = tmp_study_data['study_name']
    study_padmet = tmp_study_data['study_padmet']
    verbose = tmp_study_data['verbose']
    orthology_based_path = tmp_study_data['orthology_based_path']
    database_path = tmp_study_data['database_path']
    padmet_from_networks_path = tmp_study_data['padmet_from_networks_path']
    sbml_from_networks_path = tmp_study_data['sbml_from_networks_path']

    padmet_output = "{0}/{1}.padmet".format(padmet_from_networks_path, study_name)
    sbml_output = "{0}/{1}.sbml".format(sbml_from_networks_path, study_name)
    if os.path.exists(padmet_output):
        if verbose:
            print("%s already exist, skip" %os.path.basename(padmet_output))
            return
    else:
        ortho_sbml_folder = "{0}/{1}".format(orthology_based_path, study_name)
        source_tool = "ORTHOFINDER"
        source_category = "ORTHOLOGY"
        if verbose:
            print("Creating %s" %os.path.basename(padmet_output))
        if os.path.exists(study_padmet):
            if verbose:
                print("\tStarting from %s" %os.path.basename(study_padmet))
            padmet_path = study_padmet
            if os.path.exists(ortho_sbml_folder):
                sbml_to_padmet.sbml_to_padmetSpec(ortho_sbml_folder, padmet_path, padmetRef_file=database_path, output=padmet_output, source_tool=source_tool, source_category=source_category, verbose=verbose)
            else:
                if verbose:
                    print("\tNo orthology folder.")
                    print(("\tMove {0} in {1}".format(study_name, padmet_output)))
                copyfile(padmet_path, padmet_output)
                return
        if os.path.exists(ortho_sbml_folder) and next(os.walk(ortho_sbml_folder))[2]:
            sbml_to_padmet.sbml_to_padmetSpec(ortho_sbml_folder, padmet_output, padmetRef_file=database_path, source_tool=source_tool, source_category=source_category, verbose=verbose)
            if not os.path.isfile(sbml_output):
                if os.path.isfile(padmet_output):
                    if verbose:
                        print("Creating sbml from padmet for %s" %study_name)
                    sbmlGenerator.padmet_to_sbml(padmet=padmet_output, output=sbml_output, verbose=verbose)
                else:
                    if verbose:
                        print("\tNo padmet file to create sbml for %s'" %study_name)
                    return

            else:
                if verbose:
                    print("\t%s's sbml alreayd exists" %study_name)
                return

        else:
            if verbose:
                print("\t%s's folder is empty" %study_name)
            return
