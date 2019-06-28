#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
usage:
    aucome analysis --run=ID [--cpu=INT] [-v]

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
import pandas as pa
import re
import subprocess
import sys
import time

from aucome.utils import parse_config_file
from Bio import SeqIO
from multiprocessing import Pool


def command_help():
    print(docopt.docopt(__doc__))


def analysis_parse_args(command_args):
    args = docopt.docopt(__doc__, argv=command_args)
    run_id = args['--run']
    verbose = args['-v']

    if args["--cpu"]:
        nb_cpu_to_use = int(args["--cpu"])
    else:
        nb_cpu_to_use = 1

    run_analysis(run_id, nb_cpu_to_use, verbose)


def run_analysis(run_id, nb_cpu_to_use, verbose):

    aucome_pool = Pool(nb_cpu_to_use)

    config_data = parse_config_file(run_id)

    analysis_path = config_data['analysis_path']
    analysis_group_file_path = config_data['analysis_group_file_path']

    group_data = []
    with open(analysis_group_file_path, 'r') as group_file:
        group_reader = csv.reader(group_file, delimiter='\t')
        for row in group_reader:
            tmp_data = {'padmet_utils_path': config_data['padmet_utils_path'], 'database_path': config_data['database_path'],
                        'padmet_from_networks_path': config_data['padmet_from_networks_path'], 'analysis_path': analysis_path,
                        'studied_organisms_path': config_data['studied_organisms_path'], 'verbose': verbose, 'group_name': row[0],
                        'groups': row[1:]}

            group_data.append(tmp_data)

    aucome_pool.map(analysis_on_group, group_data)

    aucome_pool.close()
    aucome_pool.join()

def analysis_on_group(group_data):
    group_name = group_data['group_name']
    groups = group_data['groups']
    studied_organisms_path = group_data['studied_organisms_path']
    padmet_utils_path = group_data['padmet_utils_path']
    database_path = group_data['database_path']
    padmet_from_networks_path = group_data['padmet_from_networks_path']
    analysis_path = group_data['analysis_path']
    verbose = group_data['verbose']
    group_name = group_data['group_name']

    group_analysis_path = analysis_path + '/' + group_name

    if not os.path.isdir(group_analysis_path):
        if len(groups) == 1:
            sys.exit('A group must contain more than one member.')

        for species in groups:
            if species not in os.listdir(studied_organisms_path):
                sys.exit(species + ' not a valid species name from :' + '\n'.join(os.listdir(studied_organisms_path)))

        cmds = ["python3",  padmet_utils_path + "/padmet_utils/exploration/compare_padmet.py", "--padmet", padmet_from_networks_path,
                "--output", group_analysis_path, "--padmetRef", database_path]

        if verbose:
            cmds.append('-v')

        subprocess.call(cmds)

        for csv_file in os.listdir(group_analysis_path):
            tmp_df = pa.read_csv(group_analysis_path + '/' + csv_file, sep='\t')
            new_columns = [column for column in tmp_df.columns
                            if 'reaction' == column or any([species in column for species in groups])]
            tmp_df = tmp_df[new_columns]
            tmp_df.to_csv(group_analysis_path + '/' + csv_file, sep='\t', index=False)

        cmds = ["python3",  padmet_utils_path + "/padmet_utils/exploration/dendrogram_reactions_distance.py", "--reactions", group_analysis_path + '/reactions.csv',
                "--output", group_analysis_path + '/dendrogram_output', "--padmetRef", database_path]

        if verbose:
            cmds.append('-v')

        subprocess.call(cmds)

    else:
        print(group_analysis_path + ' already exists. Delete it if you want to relaunch the analysis.')