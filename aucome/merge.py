#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
usage:
    aucome merge --run=ID [--cpu=INT] [-v] [--vv]

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
import time

from shutil import copyfile
from padmet.utils.connection import sbml_to_padmet, sbmlGenerator, padmet_to_padmet
from padmet.classes import PadmetRef, PadmetSpec
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
    if verbose:
        print('--- Running merge step ---')
    merge_start_time = time.time()
    aucome_pool = Pool(nb_cpu_to_use)

    config_data = parse_config_file(run_id)

    padmet_from_annotation_path = config_data['padmet_from_annotation_path']
    padmet_from_networks_path = config_data['padmet_from_networks_path']
    sbml_from_networks_path = config_data['sbml_from_networks_path']
    database_path = config_data['database_path']

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
                            'sbml_from_networks_path': sbml_from_networks_path, 'database_path': database_path,
                            'verbose': verbose, 'veryverbose': veryverbose}
        study_draft_data.append(tmp_study_data)
    aucome_pool.map(create_output, study_draft_data)

    aucome_pool.close()
    aucome_pool.join()

    padmet_to_padmet.padmet_to_padmet(padmet_from_networks_path, networks_path + '/panmetabolism.padmet', verbose=veryverbose)
    sbmlGenerator.padmet_to_sbml(padmet=networks_path + '/panmetabolism.padmet', output=networks_path + '/panmetabolism.sbml', verbose=veryverbose)

    merge_end_time = (time.time() - merge_start_time)
    integer_part, decimal_part = str(merge_end_time).split('.')
    merge_time = ".".join([integer_part, decimal_part[:3]])

    if verbose:
        print("--- merge step done in: %ss ---" %merge_time)


def add_spontaneous_reactions(padmet_path, padmet_ref_path, output_padmet_path, only_complete_pathways=True):
    number_spontaneous_reactions = 0

    padmetSpec = PadmetSpec(padmet_path)
    padmetRef = PadmetRef(padmet_ref_path)

    all_spontaneous_rxns = set([node.id for node in list(padmetRef.dicOfNode.values()) if node.type == "reaction" and "SPONTANEOUS" in node.misc])

    for spontaneous_rxn_id in all_spontaneous_rxns:
        in_pwys = set([rlt.id_out for rlt in padmetRef.dicOfRelationIn.get(spontaneous_rxn_id,None) if rlt.type == "is_in_pathway"])
        for pwy_id in in_pwys:
            if pwy_id in padmetSpec.dicOfNode.keys():
                padmet_ref_in_rxns = set([rlt.id_in for rlt in padmetRef.dicOfRelationOut.get(pwy_id,[]) if rlt.type == "is_in_pathway"])
                padmet_spec_in_rxns = set([rlt.id_in for rlt in padmetSpec.dicOfRelationOut.get(pwy_id,[]) if rlt.type == "is_in_pathway"])

                if only_complete_pathways:
                    difference_rxns = padmet_ref_in_rxns.difference(padmet_spec_in_rxns)

                    if difference_rxns != set():
                        if difference_rxns.issubset(all_spontaneous_rxns):
                            for difference_rxn in difference_rxns:
                                if difference_rxn not in set([node.id for node in list(padmetSpec.dicOfNode.values()) if node.type == "reaction"]):
                                    padmetSpec.copyNode(padmetRef, difference_rxn)
                                    number_spontaneous_reactions += 1
                else:
                    if spontaneous_rxn_id not in set([node.id for node in list(padmetSpec.dicOfNode.values()) if node.type == "reaction"]):
                        padmetSpec.copyNode(padmetRef, spontaneous_rxn_id)
                        number_spontaneous_reactions += 1

    padmetSpec.generateFile(output_padmet_path)

    print('Add {0} spontaneous reactions to {1}'.format(number_spontaneous_reactions, output_padmet_path))


def create_output(tmp_study_data):
    padmet_path = tmp_study_data['padmet_path']
    study_padmet = tmp_study_data['study_padmet'].replace('.padmet', '').replace('output_pathwaytools_', '')
    verbose = tmp_study_data['verbose']
    veryverbose = tmp_study_data['veryverbose']
    padmet_from_networks_path = tmp_study_data['padmet_from_networks_path']
    sbml_from_networks_path = tmp_study_data['sbml_from_networks_path']
    padmet_ref_path = tmp_study_data['database_path']

    if not os.path.exists(padmet_from_networks_path + '/' + study_padmet + '.padmet'):
        if verbose:
            print('Create ' + study_padmet +' from ' + padmet_path + ' to ' + padmet_from_networks_path)
        add_spontaneous_reactions(padmet_path, padmet_ref_path, padmet_from_networks_path + '/' + study_padmet + '.padmet')
    else:
        print('There is already a padmet for ' + study_padmet + ' ' + padmet_from_networks_path + '.')

    if not os.path.exists(sbml_from_networks_path + '/' + study_padmet + '.sbml'):
        sbmlGenerator.padmet_to_sbml(padmet=padmet_path, output=sbml_from_networks_path + '/' + study_padmet + '.sbml', verbose=veryverbose)
    else:
        print('There is already a sbml for ' + study_padmet + ' in ' + sbml_from_networks_path + '.')
