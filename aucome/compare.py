#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
usage:
    aucome compare --run=ID [--cpu=INT] [-v]

options:
    --run=ID    Pathname to the comparison workspace.
    --cpu=INT     Number of cpu to use for the multiprocessing (if none use 1 cpu). [default: 1]
    -v     Verbose.

"""

import csv
import docopt
import matplotlib.pyplot as plt
import os
import pandas as pa
import time

from aucome.utils import parse_config_file
from padmet.classes.padmetRef import PadmetRef
from padmet.utils.exploration import compare_padmet, dendrogram_reactions_distance
from supervenn import supervenn


def command_help():
    """
    Show help for analysis command.
    """
    print(docopt.docopt(__doc__))


def compare_parse_args(command_args):
    """ Parse args from __main__.py

    Args:
        command_args (list): command arguments
    """
    args = docopt.docopt(__doc__, argv=command_args)
    run_id = args['--run']
    verbose = args['-v']

    if args["--cpu"]:
        nb_cpu_to_use = int(args["--cpu"])
    else:
        nb_cpu_to_use = 1

    run_compare(run_id, nb_cpu_to_use, verbose)


def run_compare(run_id, nb_cpu_to_use, verbose):
    """Compare the group specified by the user.

    Args:
        run_id (str): ID of the run
        nb_cpu_to_use (int): number of CPU for multiprocessing
        verbose (boolean): verbose
    """
    if verbose:
        print('--- Running compare step ---')
    compare_start_time = time.time()
    config_data = parse_config_file(run_id)

    analysis_path = config_data['analysis_path']
    analysis_group_file_path = config_data['analysis_group_file_path']
    compare_output_path = analysis_path + '/compare_group'

    database_path = config_data['database_path']
    padmet_from_networks_path = config_data['padmet_from_networks_path']

    # Create a dictionary containing the group name and the species inside the group.
    group_data = {}
    padmets = []
    with open(analysis_group_file_path, 'r') as group_file:
        group_reader = csv.reader(group_file, delimiter='\t')
        cluster_reactions = {}
        for row in group_reader:
            group_name = row[0]
            groups = [species for species in row[1:] if species != '']
            group_data[group_name] = groups
            if group_name != 'all':
                padmets.extend([padmet_from_networks_path + '/' + species + '.padmet' for species in groups])

    padmets = list(set(padmets))

    if not os.path.isdir(compare_output_path):
        os.mkdir(compare_output_path)

    padmetref = PadmetRef(database_path)
    # Create the reactions.tsv file needed to create dendrogram.
    padmet_path = ','.join(padmets)
    compare_padmet.compare_padmet(padmet_path=padmet_path, output=compare_output_path, padmetRef=padmetref, verbose=verbose)

    # Read the reactions.tsv file and remove the column unused.
    reactions_file = compare_output_path + '/' + 'reactions.tsv'
    reactions_dataframe = pa.read_csv(reactions_file, sep='\t')
    columns = [column for column in reactions_dataframe.columns if '(sep=;)' not in column and '_formula' not in column]
    reactions_dataframe = reactions_dataframe[columns].copy()
    reactions_dataframe.set_index('reaction', inplace=True)

    # For each group, extract the reactions present in its species to create supervenn sets.
    supervenn_sets = []
    supervenn_labels = []
    for group_name in group_data:
        if group_name != 'all':
            groups = group_data[group_name]
            reactions_temp = []
            for species in groups:
                species_reactions_dataframe = reactions_dataframe[reactions_dataframe[species] == 1]
                reactions_temp.extend(species_reactions_dataframe.index.tolist())
            supervenn_sets.append(set(reactions_temp))
            supervenn_labels.append(group_name)
            cluster_reactions[group_name] = set(reactions_temp)

    supervenn(supervenn_sets, supervenn_labels, chunks_ordering='occurrence', sets_ordering='minimize gaps')
    plt.savefig(compare_output_path + '/compare_group.png', bbox_inches='tight')
    plt.clf()

    dendrogram_reactions_distance.reaction_figure_creation(reactions_file, os.path.join(compare_output_path, "dendrogram_output"), padmetRef_file=database_path, verbose=verbose)

    compare_end_time = (time.time() - compare_start_time)
    integer_part, decimal_part = str(compare_end_time).split('.')
    compare_time = ".".join([integer_part, decimal_part[:3]])

    if verbose:
        print("--- compare step done in: %ss ---" %compare_time)
