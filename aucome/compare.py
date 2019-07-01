#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
usage:
    aucome compare --run=ID [--group=STR] [--cpu=INT] [-v]

options:
    --run=ID    Pathname to the comparison workspace.
    --group=STR    Group to compare.
    --cpu=INT     Number of cpu to use for the multiprocessing (if none use 1 cpu).
    -v     Verbose.

"""

import configparser
import csv
import docopt
import eventlet
import mpwt
import os
import pandas as pa
import re
import subprocess
import time

from aucome.utils import parse_config_file
from Bio import SeqIO
from multiprocessing import Pool


def command_help():
    print(docopt.docopt(__doc__))


def compare_parse_args(command_args):
    args = docopt.docopt(__doc__, argv=command_args)
    run_id = args['--run']
    group_to_compares = args['--group']
    verbose = args['-v']

    if args["--cpu"]:
        nb_cpu_to_use = int(args["--cpu"])
    else:
        nb_cpu_to_use = 1

    run_compare(run_id, nb_cpu_to_use, group_to_compares, verbose)


def run_compare(run_id, nb_cpu_to_use, group_to_compares, verbose):

    config_data = parse_config_file(run_id)

    analysis_path = config_data['analysis_path']
    analysis_group_file_path = config_data['analysis_group_file_path']
    upset_path = analysis_path + '/upset_graph'
    upset_tmp_data_path = analysis_path + '/upset_graph/tmp_data'

    padmet_utils_path = config_data['padmet_utils_path']
    database_path = config_data['database_path']
    padmet_from_networks_path = config_data['padmet_from_networks_path']

    if not os.path.isdir(upset_path):
        os.mkdir(upset_path)

        cmds = ["python3",  padmet_utils_path + "/padmet_utils/exploration/compare_padmet.py", "--padmet", padmet_from_networks_path,
                "--output", upset_path, "--padmetRef", database_path]

        if verbose:
            cmds.append('-v')

        subprocess.call(cmds)

    os.mkdir(upset_tmp_data_path)

    reactions_dataframe = pa.read_csv(upset_path + '/' + 'reactions.csv', sep='\t')
    columns = [column for column in reactions_dataframe.columns if '(sep=;)' not in column]
    columns = [column for column in columns if '_formula' not in column]
    reactions_dataframe = reactions_dataframe[columns].copy()
    reactions_dataframe.set_index('reaction', inplace=True)

    # Translate 'present'/(nan) data into a True/False absence-presence matrix.
    for column in reactions_dataframe.columns.tolist():
        reactions_dataframe[column] = [True if data == 'present' else False for data in reactions_dataframe[column]]

    group_data = {}
    with open(analysis_group_file_path, 'r') as group_file:
        group_reader = csv.reader(group_file, delimiter='\t')
        cluster_reactions = {}
        for row in group_reader:
            group_name = row[0]
            groups = row[1:]
            group_data[group_name] = groups

    intersect_species = []
    for group_name in group_data:
        if group_name != 'all':
            intersect_species.extend(group_data[group_name])

    for group_name in group_data:
        if group_name == 'all':
            groups = list(set(group_data[group_name]) - set(intersect_species))
            group_name = 'other_species'
        else:
            groups = group_data[group_name]
        reactions_temp = []
        for species in groups:
            species_reactions_dataframe = reactions_dataframe[reactions_dataframe[species] == True]
            reactions_temp.extend(species_reactions_dataframe.index.tolist())
        cluster_reactions[group_name] = set(reactions_temp)


        df = pa.DataFrame({group_name: list(cluster_reactions[group_name])})
        df.to_csv(upset_tmp_data_path+'/'+group_name+'.tsv', sep='\t', index=None, header=None)

    cmd = 'intervene upset -i  {0}/*.tsv --type list -o {1} --figtype svg'.format(upset_tmp_data_path, upset_path)
    if verbose:
        subprocess.call(cmd, shell=True)
    else:
        FNULL = open(os.devnull, 'w')
        subprocess.call(cmd, shell=True, stdout=FNULL, stderr=subprocess.STDOUT)