#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
usage:
    aucome merge --run=ID [--cpu=INT] [-v]

options:
    --run=ID    Pathname to the comparison workspace.
    --cpu=INT     Number of cpu to use for the multiprocessing (if none use 1 cpu).
    -v     Verbose.
    --vv    Very verbose.
"""

import docopt
import os
import shutil
import sys

from shutil import copyfile
from padmet.utils.connection import sbml_to_padmet, sbmlGenerator, padmet_to_padmet

from aucome.utils import parse_config_file
from multiprocessing import Pool


def command_help():
    print(docopt.docopt(__doc__))


def merge_parse_args(command_args):
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

    run_merge(run_id, nb_cpu_to_use, verbose, veryverbose)


def run_merge(run_id, nb_cpu_to_use, verbose, veryverbose=None):

    aucome_pool = Pool(nb_cpu_to_use)

    config_data = parse_config_file(run_id)

    padmet_from_annotation_path = config_data['padmet_from_annotation_path']
    padmet_from_networks_path = config_data['padmet_from_networks_path']
    sbml_from_networks_path = config_data['sbml_from_networks_path']

    structural_padmets_path = config_data['structural_padmets_path']
    orthofinder_filtered_path = config_data['orthofinder_filtered_path']
    orthofinder_padmet_path = config_data['orthofinder_padmet_path']
    padmet_from_annotation_path = config_data['padmet_from_annotation_path']
    networks_path = config_data['networks_path']

    structural_padmets = [padmet for padmet in os.listdir(structural_padmets_path) if padmet.endswith('.padmet')]
    orthofinder_filtered_padmets = [padmet for padmet in os.listdir(orthofinder_filtered_path) if padmet.endswith('.padmet')]
    orthofinder_padmets = [padmet for padmet in os.listdir(orthofinder_padmet_path) if padmet.endswith('.padmet')]
    pathway_tools_padmets = [padmet for padmet in os.listdir(padmet_from_annotation_path) if padmet.endswith('.padmet')]

    if len(structural_padmets) > 0:
        padmets = [(padmet, structural_padmets_path + '/' + padmet) for padmet in structural_padmets]
    elif len(orthofinder_filtered_padmets) > 0:
        padmets = [(padmet, orthofinder_filtered_path + '/' + padmet) for padmet in orthofinder_filtered_padmets]
    elif len(orthofinder_padmets) > 0:
        padmets = [(padmet, orthofinder_padmet_path + '/' + padmet) for padmet in orthofinder_padmets]
    elif len(pathway_tools_padmets) > 0:
        padmets = [(padmet, padmet_from_annotation_path + '/' + padmet) for padmet in pathway_tools_padmets]
    else:
        sys.exit('No padmets have been created, run reconstruction or workflow.')

    study_draft_data = []
    for study_name, padmet_path in padmets:
        tmp_study_data = {'padmet_path': padmet_path, 'study_padmet': study_name, 'padmet_from_networks_path': padmet_from_networks_path,
                            'sbml_from_networks_path': sbml_from_networks_path, 'verbose': verbose, 'veryverbose': veryverbose}
        study_draft_data.append(tmp_study_data)
    aucome_pool.map(create_output, study_draft_data)

    aucome_pool.close()
    aucome_pool.join()

    padmet_to_padmet.padmet_to_padmet(padmet_from_networks_path, networks_path + '/panmetabolism.padmet', verbose=veryverbose)
    sbmlGenerator.padmet_to_sbml(padmet=networks_path + '/panmetabolism.padmet', output=networks_path + '/panmetabolism.sbml', verbose=veryverbose)


def create_output(tmp_study_data):
    padmet_path = tmp_study_data['padmet_path']
    study_padmet = tmp_study_data['study_padmet'].replace('.padmet', '').replace('output_pathwaytools_', '')
    verbose = tmp_study_data['verbose']
    veryverbose = tmp_study_data['veryverbose']
    padmet_from_networks_path = tmp_study_data['padmet_from_networks_path']
    sbml_from_networks_path = tmp_study_data['sbml_from_networks_path']

    if not os.path.exists(padmet_from_networks_path + '/' + study_padmet + '.padmet'):
        if verbose:
            print('Move ' + study_padmet +' from ' + padmet_path + ' to ' + padmet_from_networks_path)
        shutil.copyfile(padmet_path, padmet_from_networks_path + '/' + study_padmet + '.padmet')
    else:
        print('There is already a padmet for ' + study_padmet + ' ' + padmet_from_networks_path + '.')

    if not os.path.exists(sbml_from_networks_path + '/' + study_padmet + '.sbml'):
        sbmlGenerator.padmet_to_sbml(padmet=padmet_path, output=sbml_from_networks_path + '/' + study_padmet + '.sbml', verbose=veryverbose)
    else:
        print('There is already a sbml for ' + study_padmet + ' in ' + sbml_from_networks_path + '.')
