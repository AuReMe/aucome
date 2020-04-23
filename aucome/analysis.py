#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
usage:
    aucome analysis --run=ID [--cpu=INT] [--pvclust] [-v]

options:
    --run=ID    Pathname to the comparison workspace.
    --cpu=INT     Number of cpu to use for the multiprocessing (if none use 1 cpu). [default: 1]
    --pvclust     Run also pvclust to create another dendrogram of reaction.
    -v     Verbose.

"""

import csv
import docopt
import os
import sys
import time

from multiprocessing import Pool

from padmet.utils.connection import sbmlGenerator, padmet_to_padmet
from padmet.utils.exploration import compare_padmet, dendrogram_reactions_distance
from padmet.classes import PadmetRef

from aucome.utils import parse_config_file


def command_help():
    """
    Show help for analysis command.
    """
    print(docopt.docopt(__doc__))


def analysis_parse_args(command_args):
    """ Parse args from __main__.py

    Args:
        command_args (list): command arguments
    """
    args = docopt.docopt(__doc__, argv=command_args)
    run_id = args['--run']
    verbose = args['-v']
    nb_cpu_to_use = int(args["--cpu"])
    pvclust = args['--pvclust']

    run_analysis(run_id, nb_cpu_to_use, pvclust, verbose)


def run_analysis(run_id, nb_cpu_to_use, pvclust, verbose):
    """Create input data for creationf of reaction dendrogram tsv reactions files.

    Args:
        run_id (str): ID of the run
        nb_cpu_to_use (int): number of CPU for multiprocessing
        pvclust (boolean): use also pvclust to create reaction dendrogram
        verbose (boolean): verbose
    """
    if verbose:
        print('--- Running analysis step ---')
    analysis_start_time = time.time()
    config_data = parse_config_file(run_id)

    analysis_group_file_path = config_data['analysis_group_file_path']

    # Create list of dictionaries containing input data for multiprocessing.
    # As we have to give one argument after the function to pool().
    # Create one dictionary by group in the group_template.tsv file.
    group_data = []
    with open(analysis_group_file_path, 'r') as group_file:
        group_reader = csv.reader(group_file, delimiter='\t')
        for row in group_reader:
            group_name = row[0]
            groups = [org_name for org_name in row[1:] if org_name]
            analysis_on_group(group_name, groups, config_data, pvclust, nb_cpu_to_use, verbose)

    analysis_end_time = (time.time() - analysis_start_time)
    integer_part, decimal_part = str(analysis_end_time).split('.')
    analysis_time = ".".join([integer_part, decimal_part[:3]])

    if verbose:
        print("--- analysis step done in: %ss ---" %analysis_time)


def analysis_on_group(group_name, groups, config_data, pvclust, nb_cpu_to_use, verbose):
    """Create reaction dendrogram and extract specific reactions using metabolic networks.

    Args:
        group_name (str): Name of the group from group_template.tsv.
        groups (list): All the species inside the group.
        config_data (dict): Dictionary with all configuration paths.
        pvclust (boolean): use also pvclust to create reaction dendrogram
        nb_cpu_to_use (int): number of CPU for multiprocessing
        verbose (bool): Verbose.
    """

    database_path = config_data['database_path']
    padmetRef = PadmetRef(database_path)
    padmet_from_networks_path = config_data['padmet_from_networks_path']
    analysis_path = config_data['analysis_path']

    all_padmet_path = [os.path.join(padmet_from_networks_path,name+".padmet") for name in groups ]
    group_analysis_path = analysis_path + '/' + group_name

    if not os.path.isdir(group_analysis_path):
        if len(groups) == 1:
            sys.exit('A group must contain more than one member.')

        for padmet_path in all_padmet_path:
            if not os.path.exists(padmet_path):
                org_name = os.path.splitext(os.path.basename(padmet_path))[0]
                sys.exit("Padmet file of organism %s from group %s not found in %s" %(org_name, group_name, padmet_from_networks_path))

        # Compare the padmet to create the reactions.tsv file needed to create the reaction dendrogram.
        compare_padmet.compare_padmet(padmet_path=",".join(all_padmet_path), output=group_analysis_path, padmetRef=padmetRef, verbose=verbose, number_cpu=nb_cpu_to_use)
        padmet_to_padmet.padmet_to_padmet(",".join(all_padmet_path), group_analysis_path + '/' + group_name + '_panmetabolism.padmet')
        sbmlGenerator.padmet_to_sbml(padmet=group_analysis_path + '/' + group_name + '_panmetabolism.padmet', output=group_analysis_path + '/' + group_name + '_panmetabolism.sbml', verbose=verbose)

        dendrogram_reactions_distance.reaction_figure_creation(reaction_file=group_analysis_path + '/reactions.tsv', output_folder=group_analysis_path + '/dendrogram_output',
                                                                padmetRef_file=database_path, pvclust=pvclust, verbose=verbose)

    else:
        print(group_analysis_path + ' already exists. Delete it if you want to relaunch the analysis.')