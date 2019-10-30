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
import os
import pandas as pa
import subprocess

from padmet.utils.exploration import compare_padmet, dendrogram_reactions_distance

from aucome.utils import parse_config_file


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
    """Compare the gorup specified by the user.

    Args:
        run_id (str): ID of the run
        nb_cpu_to_use (int): number of CPU for multiprocessing
        verbose (boolean): verbose
    """
    config_data = parse_config_file(run_id)

    analysis_path = config_data['analysis_path']
    analysis_group_file_path = config_data['analysis_group_file_path']
    upset_path = analysis_path + '/upset_graph'
    upset_tmp_data_path = upset_path + '/tmp_data'
    upset_tmp_reaction_path = upset_tmp_data_path + '/tmp'

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

    if not os.path.isdir(upset_path):
        os.mkdir(upset_path)
    if not os.path.isdir(upset_tmp_data_path):
        os.mkdir(upset_tmp_data_path)
        if not os.path.isdir(upset_tmp_reaction_path):
            os.mkdir(upset_tmp_reaction_path)

        # Create the reactions.csv file needed to create dendrogram.
        padmet_path = ','.join(padmets)
        compare_padmet.compare_padmet(padmet_path, upset_tmp_reaction_path, database_path, verbose)

    # Read the reactions.csv file and remove the column unused.
    reactions_file = upset_tmp_reaction_path + '/' + 'reactions.csv'
    reactions_dataframe = pa.read_csv(reactions_file, sep='\t')
    columns = [column for column in reactions_dataframe.columns if '(sep=;)' not in column]
    columns = [column for column in columns if '_formula' not in column]
    reactions_dataframe = reactions_dataframe[columns].copy()
    reactions_dataframe.set_index('reaction', inplace=True)

    # Translate 'present'/(nan) data into a True/False absence-presence matrix.
    for column in reactions_dataframe.columns.tolist():
        reactions_dataframe[column] = [True if data == 'present' else False for data in reactions_dataframe[column]]

    # For each group, extract the reactions present in its species.
    # Then create a tsv file containing these reactions..
    for group_name in group_data:
        if group_name != 'all':
            groups = group_data[group_name]
            reactions_temp = []
            for species in groups:
                species_reactions_dataframe = reactions_dataframe[reactions_dataframe[species] == True]
                reactions_temp.extend(species_reactions_dataframe.index.tolist())
            cluster_reactions[group_name] = set(reactions_temp)

            df = pa.DataFrame({group_name: list(cluster_reactions[group_name])})
            df.to_csv(upset_tmp_data_path+'/'+group_name+'.tsv', sep='\t', index=None, header=None)

    # Launch Intervene to create upset graph using each group file.
    upset_data_path = [upset_tmp_data_path + '/' + tsv_file for tsv_file in os.listdir(upset_tmp_data_path) if tsv_file.endswith('.tsv')]
    cmds = ['intervene', 'upset', '-i', *upset_data_path, '--type', 'list', '-o', upset_path, '--figtype', 'svg']

    if verbose:
        subprocess.call(cmds)
    else:
        FNULL = open(os.devnull, 'w')
        subprocess.call(cmds, stdout=FNULL, stderr=subprocess.STDOUT)

    dendrogram_reactions_distance.reaction_figure_creation(reactions_file, os.path.join(upset_path, "dendrogram_output"), padmet_ref_file=database_path, verbose=verbose)

