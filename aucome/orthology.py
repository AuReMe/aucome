#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
usage:
    aucome orthology --run=ID [-S=STR] [--orthogroups] [--cpu=INT] [-v]

options:
    --run=ID    Pathname to the comparison workspace.
    --orthogroups    Use Orthogroups instead of Orthologues after Orthofinder.
    -S=STR    Sequence search program for Orthofinder [Default: diamond].
        Options: blast, mmseqs, blast_gz, diamond
    --cpu=INT     Number of cpu to use for the multiprocessing (if none use 1 cpu).
    -v     Verbose.
"""

import docopt
import os
import subprocess
import time

from padmet.utils.exploration import convert_sbml_db
from padmet.utils.connection import extract_orthofinder

from aucome.utils import parse_config_file
from multiprocessing import Pool


def command_help():
    print(docopt.docopt(__doc__))


def orthology_parse_args(command_args):
    args = docopt.docopt(__doc__, argv=command_args)
    run_id = args['--run']
    orthogroups = args['--orthogroups']
    sequence_search_prg = args['-S']
    verbose = args['-v']

    if args["--cpu"]:
        nb_cpu_to_use = int(args["--cpu"])
    else:
        nb_cpu_to_use = 1

    run_orthology(run_id, orthogroups, sequence_search_prg, nb_cpu_to_use, verbose)

def run_orthology(run_id, orthogroups, sequence_search_prg, nb_cpu_to_use, verbose):
    aucome_pool = Pool(nb_cpu_to_use)

    config_data = parse_config_file(run_id)

    orthofinder_wd_path = config_data['orthofinder_wd_path']
    orthofinder_bin_path = config_data['orthofinder_bin_path']
    orthology_based_path = config_data['orthology_based_path']
    studied_organisms_path = config_data['studied_organisms_path']
    model_organisms_path = config_data['model_organisms_path']
    mnx_cpd_path = config_data['mnx_cpd_path']
    mnx_rxn_path = config_data['mnx_rxn_path']

    all_study_name = set(next(os.walk(studied_organisms_path))[1])
    all_model_name = set(next(os.walk(model_organisms_path))[1])

    all_study_faa = dict([(study_name, "{0}/{1}/{1}.faa".format(studied_organisms_path, study_name))
                          if os.path.isfile("{0}/{1}/{1}.faa".format(studied_organisms_path, study_name))
                          else (study_name, '')
                          for study_name in all_study_name])

    all_model_faa = dict([(model_name, "{0}/{1}/{1}.faa".format(model_organisms_path, model_name))
                          if os.path.isfile("{0}/{1}/{1}.faa".format(model_organisms_path, model_name))
                          else (model_name, '')
                          for model_name in all_model_name])

    #check if Orthofinder already run, if yes, get the last workdir
    try:
        if orthogroups:
            orthodata_path = max(["%s/%s" %(x[0], 'Orthogroups/Orthogroups.tsv') for x in os.walk(orthofinder_wd_path) if 'Orthogroups' in x[1]])
        else:
            orthodata_path = max(["%s/%s" %(x[0], 'Orthologues') for x in os.walk(orthofinder_wd_path) if 'Orthologues' in x[1]])
    except ValueError:
        if verbose:
            print("Enable to find file Orthogroups.csv in {0}, need to run Orthofinder...".format(orthofinder_wd_path))
        for name, faa_path in list(all_study_faa.items()):
            if not os.path.isfile("{0}/{1}.faa".format(orthofinder_wd_path, name)):
                if verbose:
                    print("Copying {0}'s faa to {1}".format(name, orthofinder_wd_path))
                cmds = ["cp", faa_path, orthofinder_wd_path]
                subprocess.call(cmds)
        for name, faa_path in list(all_model_faa.items()):
            if not os.path.isfile("{0}/{1}.faa".format(orthofinder_wd_path, name)):
                if verbose:
                    print("Copying {0}'s faa to {1}".format(name, orthofinder_wd_path))
                cmds = ["cp", faa_path, orthofinder_wd_path]
                subprocess.call(cmds)

        if verbose:
            print("Running Orthofinder on %s cpu" %nb_cpu_to_use)

        chronoDepart = time.time()
        cmds = [orthofinder_bin_path, "-f", orthofinder_wd_path, "-t", str(nb_cpu_to_use),
                "-S", sequence_search_prg]
        subprocess.call(cmds)
        chrono = (time.time() - chronoDepart)
        partie_entiere, partie_decimale = str(chrono).split('.')
        chrono = ".".join([partie_entiere, partie_decimale[:3]])
        if verbose:
            print("Orthofinder done in: %ss" %chrono)
        if orthogroups:
            orthodata_path = max(["%s/%s" %(x[0], 'Orthogroups/Orthogroups.tsv') for x in os.walk(orthofinder_wd_path) if 'Orthogroups' in x[1]])
        else:
            orthodata_path = max(["%s/%s" %(x[0], 'Orthologues') for x in os.walk(orthofinder_wd_path) if 'Orthologues' in x[1]])
    if verbose:
        print("Parsing Orthofinder output %s" %orthodata_path)

    if verbose:
        print("Start sbml creation...")
    all_dict_data = []
    for study_name in all_study_name:
        dict_data = {'sbml': run_id, 'orthodata_path': orthodata_path, 'study_name': study_name,
                    'verbose': verbose, 'orthogroups': orthogroups,
                    'output': orthology_based_path + '/' + study_name}
        all_dict_data.append(dict_data)

    chronoDepart = time.time()
    aucome_pool.map(orthogroup_to_sbml, all_dict_data)
    chrono = (time.time() - chronoDepart)
    integer_part, decimal_part = str(chrono).split('.')
    chrono = ".".join([integer_part, decimal_part[:3]])
    if verbose:
        print("Orthofinder output parsed in: %ss" %chrono)
    """
    #check database, mapping to metacyc ???
    data_convert_sbml_db = []
    for dict_data in all_dict_data:
        tmp_dict_data = {'sbml': orthology_based_path + '/' + study_name, 
                         'mnx_rxn_path': mnx_rxn_path, 'mnx_cpd_path': mnx_cpd_path, 'verbose': verbose}
        data_convert_sbml_db.append(tmp_dict_data)
        
    aucome_pool.map(_convert_sbml_db, data_convert_sbml_db)

    aucome_pool.close()
    aucome_pool.join()
    """

def _convert_sbml_db(data_convert_sbml_db):
    
    sbml_file = data_convert_sbml_db['sbml']
    verbose = data_convert_sbml_db['verbose']
    mnx_rxn_path = data_convert_sbml_db['mnx_rxn_path']
    mnx_cpd_path = data_convert_sbml_db['mnx_cpd_path']

    if os.path.isfile(sbml_file):
        dict_file = "{0}_dict.csv".format(os.path.splitext(sbml_file)[0])
        if not os.path.exists(dict_file):
            db_ref = convert_sbml_db.check_sbml_db(sbml_file, "reaction", mnx_reac_file=mnx_rxn_path, verbose=True)[0]
            if verbose:
                print("%s: %s" %(os.path.basename(sbml_file), db_ref))
            if db_ref.lower() != "metacyc":
                if verbose:
                    print("Creating id mapping file: %s" %dict_file)
                convert_sbml_db.map_sbml(sbml_file, "reaction", "metacyc", dict_file, verbose=verbose, mnx_rxn_path=mnx_rxn_path, mnx_chem_path=mnx_cpd_path)


def orthogroup_to_sbml(dict_data):
    """
    dict_orthogroup: global var
    """
    #dict_data = {'study_name':'', 'o_compare_name': '', sbml_template':'', 'output':''}
    sbml = dict_data['sbml']
    orthodata_path = dict_data['orthodata_path']
    study_name = dict_data['study_name']
    output = dict_data['output']
    verbose = dict_data['verbose']
    orthogroups = dict_data['orthogroups']

    all_model_sbml = extract_orthofinder.get_sbml_files(sbml, workflow="aucome", verbose=verbose)
    if orthogroups:
        extract_orthofinder.orthogroups_to_sbml(orthodata_path, all_model_sbml, output, study_name, verbose)
    else:
        extract_orthofinder.orthologue_to_sbml(orthodata_path, all_model_sbml, output, study_name, verbose)
