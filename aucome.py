#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Apr 26 16:32:23 2018

@author: Meziane AITE, meziane.aite@inria.fr

1./ get all fasta in studied_orgnisms/ in dict_faa_paths['study']
ex dict_faa_paths['study'] = [studied_organisms_path/study_1/faa_study_name, ...]
get all fasta in model_data/
ex dict_faa_paths['model'] = [model_data_path/model_1/faa_model_name, ...]

usage:
    aucome.py --setWorkingFolder=DIR
    aucome.py --init=ID [-v]
    aucome.py --run=DIR [-c] [-o] [-S=STR] [-p] [-d] [--cpu=INT] [-v] [--log=FILE]
    aucome.py -R
    aucome.py --version
    aucome.py --installPWT=PWT_path [--ptools=ptools_path]
    aucome.py --uninstallPWT

options:
    -h --help     Show help.
    -R     Open access from container.
    --run=ID    pathname to the comparison workspace
    -c    Check inputs validity
    -o    Run Orthofinder
    -S=STR    Sequence search program for Orthofinder [Default: diamond].
          Options: blast, mmseqs, blast_gz, diamond
    -p    Run Pathway-Tools
    -d    Run Orthofinder, Pathway and merge all networks
    --cpu=INT     Number of cpu to use for the multiprocessing (if none use 1 cpu).
    -v     Verbose.

"""
__version__ = "0.4"

import configparser
import csv
import docopt
import eventlet
import itertools
import libsbml
import mpwt
import os
import re
import subprocess
import time

from Bio import SeqIO
from Bio.SeqFeature import FeatureLocation
from multiprocessing import Pool, cpu_count
from padmet.utils import sbmlPlugin as sp
from padmet.utils import gbr
from padmet.classes import PadmetSpec, PadmetRef

# Monkey patch before requests import to avoid RecursionError.
# See https://github.com/gevent/gevent/issues/1016
#eventlet.monkey_patch()
import requests


def main():
    args = docopt.docopt(__doc__)

    global all_run_folder, database_path, studied_organisms_path, model_data_path, orthology_based_path, annotation_based_path,\
    sbml_from_annotation_path, networks_path, orthofinder_bin_path, mnx_rxn_path, mnx_cpd_path,\
    sbml_study_prefix, all_study_name, all_study_padmet, all_mode_name, padmet_utils_path, release_on_gitlab, verbose,\
    all_study_gbk, model_organisms_path, all_model_gbk, padmet_from_annotation_path, study_from_annot_prefix,\
    all_study_pgdb

    release_on_gitlab = "https://gitlab.inria.fr/DYLISS/compare_metabo/raw/master/release.txt"

    # Variable of the working directory modified by setWorkingFolder arguments (modify_working_folder function).
    all_run_folder = "/shared"

    """
    args = {"--run":"test", "-v":True}
    all_run_folder = "/home/maite/Forge/docker/comparison_workspace/workdir"
    database_path = "/home/maite/Forge/docker/comparison_workspace/folder_git/compare_metabo/database/BIOCYC/METACYC/22.0_enhanced/metacyc_22.0_enhanced.padmet"
    padmet_utils_path = "/home/maite/Aureme/padmet-utils"
    mnx_rxn_path = "/home/maite/Forge/docker/comparison_workspace/folder_git/compare_metabo/database/MNX/reac_xref.tsv"
    mnx_cpd_path = "/home/maite/Forge/docker/comparison_workspace/folder_git/compare_metabo/database/MNX/chem_xref.tsv"
    """

    if args["--version"]:
        online_version = get_version()
        current_version = __version__
        if online_version:
            print("You are using the version %s, the latest is %s" %(current_version, online_version))
        else:
            print('No internet connection. Skip checking AuCoMe version.')
        return

    if args["--setWorkingFolder"]:
        modify_working_folder(args["--setWorkingFolder"])
        return

    #always_check_version
    online_version = get_version()
    current_version = __version__
    if online_version:
        if online_version != current_version:
            print("/!\ WARNING, your AuCoMe is not up-to-date. You are using the version %s, the latest is %s" %(current_version, online_version))
            print("Check the Changelog here %s" %release_on_gitlab.replace("/raw/","/blob/"))
            print("To update AuCoMe:")
            print("\tRemove your compare-img and the container created from this image:")
            print("\t\t$sudo docker rmi -f docker.io/dyliss/compare-img")
            print("\tCreate a new container with the new image:")
            print("\t\t$sudo docker run -ti -v /PATH/TO/COMPARE_WORKSPACE:/shared --name=compare docker.io/dyliss/compare-img bash")
            print("\tObviously change /PATH/TO/COMPARE_WORKSPACE to the real path of you Compare workspace")

    #add permission to all folder in all_run_folder, usefull because all cmd exec from container are root based
    if args['-R']:
        chmod_cmds = ["chmod", "-R", "777", all_run_folder]
        subprocess.call(chmod_cmds)
        return

    if args['--installPWT']:
        installing_pwt(args['--installPWT'], args['--ptools'])
        return

    if args['--uninstallPWT']:
        uninstalling_pwt()
        return

    if args["--run"]:
        run_id = args["--run"]
    else:
        run_id = args["--init"]
        
    config_file_path = "{0}/{1}/config.txt".format(all_run_folder, run_id)

    if args["--init"]:
        create_run(run_id)
        chmod_cmds = ["chmod", "-R", "777", all_run_folder]
        subprocess.call(chmod_cmds)
        return

    config = configparser.ConfigParser()
    config.read(config_file_path)
    #config.read('/home/maite/Forge/docker/comparison_workspace/template/config.txt')
    #DATABASE_PATHS
    database_path = config.get('DATABASE_PATHS','database_ref_path')
    mnx_rxn_path = config.get('DATABASE_PATHS','mnx_rxn_path')
    mnx_cpd_path = config.get('DATABASE_PATHS','mnx_cpd_path')
    #PATHS_IN_RUN
    studied_organisms_path = "{0}/{1}".format(all_run_folder, config.get('PATHS_IN_RUN','studied_organisms_path'))
    model_organisms_path = "{0}/{1}".format(all_run_folder, config.get('PATHS_IN_RUN','model_organisms_path'))
    orthology_based_path = "{0}/{1}".format(all_run_folder, config.get('PATHS_IN_RUN','orthology_based_path'))
    orthofinder_wd_path = "{0}/{1}".format(all_run_folder, config.get('PATHS_IN_RUN','orthofinder_wd_path'))    
    annotation_based_path = "{0}/{1}".format(all_run_folder, config.get('PATHS_IN_RUN','annotation_based_path'))
    pgdb_from_annotation_path = "{0}/{1}".format(all_run_folder, config.get('PATHS_IN_RUN','pgdb_from_annotation_path'))
    padmet_from_annotation_path = "{0}/{1}".format(all_run_folder, config.get('PATHS_IN_RUN','padmet_from_annotation_path'))
    sbml_from_annotation_path = "{0}/{1}".format(all_run_folder, config.get('PATHS_IN_RUN','sbml_from_annotation_path'))
    networks_path = "{0}/{1}".format(all_run_folder, config.get('PATHS_IN_RUN','networks_path'))
    #TOOL_PATHS
    orthofinder_bin_path = config.get('TOOL_PATHS','orthofinder_bin_path')
    padmet_utils_path = config.get('TOOL_PATHS', 'padmet_utils_path')
    #VAR
    study_from_annot_prefix = config.get('VAR','study_from_annot_prefix')

    if args["-v"]:
        print("Verbose Mode:")
        verbose = '-v'
    else:
        verbose = ''

    #create dict for ortho data
    all_study_name = set(next(os.walk(studied_organisms_path))[1])
    all_model_name = set(next(os.walk(model_organisms_path))[1])
    all_study_pgdb = dict([(study_name, "{0}/{1}".format(pgdb_from_annotation_path, study_name))
                          if os.path.isdir("{0}/{1}".format(pgdb_from_annotation_path, study_name))
                          else (study_name, '')
                          for study_name in all_study_name])
    all_study_gbk = dict([(study_name, "{0}/{1}/{1}.gbk".format(studied_organisms_path, study_name))
                          if os.path.isfile("{0}/{1}/{1}.gbk".format(studied_organisms_path, study_name))
                          else (study_name, '')
                          for study_name in all_study_name])
    #k = folder_name in model_organisms_path, v = path to faa in this folder, faa name should be folder_name.faa
    all_model_gbk = dict([(model_name, "{0}/{1}/{1}.gbk".format(model_organisms_path, model_name))
                          if os.path.isfile("{0}/{1}/{1}.gbk".format(model_organisms_path, model_name))
                          else (model_name, '')
                          for model_name in all_model_name])
    #k = folder_name in studied_org_path, v = path to faa in this folder, faa name should be folder_name.faa
    all_study_faa = dict([(study_name, "{0}/{1}/{1}.faa".format(studied_organisms_path, study_name))
                          if os.path.isfile("{0}/{1}/{1}.faa".format(studied_organisms_path, study_name))
                          else (study_name, '')
                          for study_name in all_study_name])
    #k = folder_name in model_organisms_path, v = path to faa in this folder, faa name should be folder_name.faa
    all_model_faa = dict([(model_name, "{0}/{1}/{1}.faa".format(model_organisms_path, model_name))
                          if os.path.isfile("{0}/{1}/{1}.faa".format(model_organisms_path, model_name))
                          else (model_name, '')
                          for model_name in all_model_name])

    all_study_padmet = dict([(study_name, "{0}/{1}{2}.padmet".format(padmet_from_annotation_path, study_from_annot_prefix, study_name))
                          if os.path.isfile("{0}/{1}{2}.padmet".format(padmet_from_annotation_path, study_from_annot_prefix, study_name))
                          else (study_name, '')
                          for study_name in all_study_name])
    #sbml of study are obtained from annotation, they should be in sbml_from_annotation_path
    #k = study_name (== folder_name in studied_org_path or obtained from sbml name), v = path to sbml, sbml_study_prefi+study_name+.sbml
    all_study_sbml = dict([(study_name, "{0}/{1}{2}.sbml".format(sbml_from_annotation_path, study_from_annot_prefix, study_name))
                           if os.path.isfile("{0}/{1}{2}.sbml".format(sbml_from_annotation_path, study_from_annot_prefix, study_name))
                           else (study_name, '')
                           for study_name in all_study_name])

    #k = folder_name in model_organisms_path, v = path to sbml in this folder, sbml name should be folder_name.sbml
    all_model_sbml = dict([(model_name, "{0}/{1}/{1}.sbml".format(model_organisms_path, model_name))
                          if os.path.isfile("{0}/{1}/{1}.sbml".format(model_organisms_path, model_name))
                          else (model_name, '')
                          for model_name in all_model_name])
    #PGDB, padmet, sbml
    all_study_pgdb = dict([(study_name, "{0}/{1}".format(pgdb_from_annotation_path, study_name))
                          if os.path.isdir("{0}/{1}".format(pgdb_from_annotation_path, study_name))
                          else (study_name, '')
                          for study_name in all_study_name])

    if args["--cpu"]:
        nb_cpu_to_use = int(args["--cpu"])
    else:
        nb_cpu_to_use = 1

    aucome_pool = Pool(nb_cpu_to_use)

    if args["-c"]:
        if verbose:
            print('Checking genbank file.')
        aucome_pool.map(check_create_faa, all_study_name)

        aucome_pool.map(create_faa_model, all_model_name)

        aucome_pool.map(create_padmet_from_pgdb, all_study_name)

        aucome_pool.map(create_sbml, all_study_padmet)

    if verbose:
        print("Input summary:")
        print("* %s Studied organims:" %(len(all_study_name)))
        for study_name in all_study_name:
            print("%s:" %study_name)
            if all_study_gbk[study_name]:
                print("\tGBK: OK")
            else:
                print("\t[WARNING] No GBK found, should be in {1}/{0}/{0}.gbk".format(study_name, studied_organisms_path))
            if all_study_pgdb[study_name]:
                print("\tPGDB: OK")
            else:
                print("\t[WARNING] No PGDB found, should be in {1}/{0}".format(study_name, pgdb_from_annotation_path))
            if all_study_padmet[study_name]:
                print("\tPADMET: OK")
            else:
                print("\t[WARNING] No PADMET found, should be in {1}/{2}{0}.padmet".format(study_name, padmet_from_annotation_path, study_from_annot_prefix))
            if all_study_faa[study_name]:
                print("\tFAA: OK")
            else:
                print("\t[WARNING] No FAA found, should be in {1}/{0}/{0}.faa".format(study_name, studied_organisms_path))
            if all_study_sbml[study_name]:
                print("\tSBML: OK")
            else:
                print("\t[WARNING] No SBML found, should be in {1}/{2}{0}.sbml".format(study_name, sbml_from_annotation_path, study_from_annot_prefix))
        print("* %s models organims:" %(len(all_model_name)))
        for model_name in all_model_name:
            print("%s:" %model_name)
            if all_model_faa[model_name]:
                print("\tFAA: OK")
            else:
                print("\t[WARNING] No FAA found, should be in {1}/{0}/{0}.faa".format(model_name, model_organisms_path))
            if all_model_sbml[model_name]:
                print("\tSBML: OK")
            else:
                print("\t[WARNING] No SBML found, should be in {1}/{0}/{0}.faa".format(model_name, model_organisms_path))

    if args["-p"]:
        #check for each study if exist PGDB folder in PGDBs folder, if missing RUN ptools
        chronoDepart = time.time()

        mpwt.multiprocess_pwt(input_folder=studied_organisms_path,
                                output_folder=pgdb_from_annotation_path,
                                patho_inference=True,
                                dat_creation=True,
                                dat_extraction=True,
                                number_cpu=nb_cpu_to_use,
                                verbose=verbose)

        chrono = (time.time() - chronoDepart)
        partie_entiere, partie_decimale = str(chrono).split('.')
        chrono = ".".join([partie_entiere, partie_decimale[:3]])
        if os.listdir(pgdb_from_annotation_path) == []:
            print('Pathway-Tools inference failed!')
            return
        if verbose:
            print("Pathway-Tools done in: %ss" %chrono)

    #for each faa, check if already in ortho_based
    if args["-o"]:
        sequence_search_prg = args['-S']
        #check if Orthofinder already run, if yes, get the last workdir
        try:
            last_orthogroup_file = max(["%s/%s" %(x[0], 'Orthogroups.csv') for x in os.walk(orthofinder_wd_path) if 'Orthogroups.csv' in x[2]])
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
            cmds = [orthofinder_bin_path, "-f", orthofinder_wd_path, "-t", nb_cpu_to_use,
                    "-S", sequence_search_prg]
            subprocess.call(cmds)
            chrono = (time.time() - chronoDepart)
            partie_entiere, partie_decimale = str(chrono).split('.')
            chrono = ".".join([partie_entiere, partie_decimale[:3]])
            if verbose:
                print("Orthofinder done in: %ss" %chrono)
            last_orthogroup_file = max(["%s/%s" %(x[0], 'Orthogroups.csv') for x in os.walk(orthofinder_wd_path) if 'Orthogroups.csv' in x[2]])
        if verbose:
            print("Parsing Orthofinder output %s" %last_orthogroup_file)
        #k=orthogroup_id, v = {k = name, v = set of genes}
        dict_orthogroup = {}
        with open(last_orthogroup_file, 'r') as csvfile:
            reader = csv.DictReader(csvfile, delimiter = "\t")
            for row in reader:
                orth_id = row['']
                row.pop('')
                new_dict = dict([(name, set(genes.split(","))) for (name, genes) in list(row.items()) if genes])
                dict_orthogroup[orth_id] = new_dict

        if verbose:
            print("Start sbml creation...")
        all_dict_data = []
        for study_name in all_study_name:
            if verbose:
                print("%s:" %study_name)
            ortho_sbml_folder = "{0}/{1}".format(orthology_based_path, study_name)
            if not os.path.exists(ortho_sbml_folder):
                if verbose:
                    print("\tCreating folder %s" %ortho_sbml_folder)
                os.makedirs(ortho_sbml_folder)
            all_to_compare = all_study_name.union(all_model_name) - set([study_name])
            for to_compare_name in all_to_compare:
                output = "{0}/{1}/output_orthofinder_{1}_from_{2}.sbml".format(orthology_based_path, study_name, to_compare_name)
                try:
                    sbml_template = all_model_sbml[to_compare_name]
                except KeyError:
                    sbml_template = all_study_sbml[to_compare_name]
                if sbml_template:
                    dict_data = {'study_name': study_name, 'to_compare_name': to_compare_name, 'sbml_template': sbml_template, 'output': output, 'verbose':verbose, 'ortho': dict_orthogroup}
                    all_dict_data.append(dict_data)

        chronoDepart = time.time()
        aucome_pool.map(orthogroup_to_sbml, all_dict_data)
        chrono = (time.time() - chronoDepart)
        partie_entiere, partie_decimale = str(chrono).split('.')
        chrono = ".".join([partie_entiere, partie_decimale[:3]])
        if verbose:
            print("Orthofinder output parsed in: %ss" %chrono)
        #check database, mapping to metacyc ???
        all_sbml_from_ortho = [dict_data['output'] for dict_data in all_dict_data]
        aucome_pool.map(convert_sbml_db, all_sbml_from_ortho)

    if args["-d"]:
        aucome_pool.map(create_draft, all_study_name)

    aucome_pool.close()
    aucome_pool.join()


def check_create_faa(study_name):
    checking_genbank(study_name)
    #create Faa from gbk if no faa found
    faa_path = "{0}/{1}/{1}.faa".format(studied_organisms_path, study_name)
    gbk_file = all_study_gbk[study_name]
    if not os.path.isfile(faa_path) and gbk_file:
        if verbose:
            print("Creating faa from gbk for %s" %study_name)
        cmds = ["python3",  padmet_utils_path + "/padmet_utils/connection/gbk_to_faa.py", "--gbk", gbk_file, "--output", faa_path]
        subprocess.call(cmds)


def create_faa_model(model_name):
    faa_path = "{0}/{1}/{1}.faa".format(model_organisms_path, model_name)
    gbk_file = all_model_gbk[model_name]
    if not os.path.isfile(faa_path) and gbk_file:
        if verbose:
            print("Creating faa from gbk for %s" %model_name)
        cmds = ["python3", padmet_utils_path + "/padmet_utils/connection/gbk_to_faa.py", "--gbk", gbk_file, "--output", faa_path]
        subprocess.call(cmds)


def create_padmet_from_pgdb(study_name):
    padmet_file = "{0}/{1}{2}.padmet".format(padmet_from_annotation_path, study_from_annot_prefix, study_name)
    pgdb_folder = all_study_pgdb[study_name]
    if not os.path.isfile(padmet_file) and pgdb_folder:
        if verbose:
            print("Creating padmet from pgdb for %s" %study_name)
        cmds = ["python3",  padmet_utils_path + "/padmet_utils/connection/pgdb_to_padmet.py", "--pgdb", pgdb_folder, "--output", padmet_file,
                "--padmetRef", database_path, "--source=genome", "--extract-gene", "--no-orphan", verbose]
        subprocess.call(cmds)


def create_sbml(study_name):
    sbml_file = "{0}/{1}{2}.sbml".format(sbml_from_annotation_path, study_from_annot_prefix, study_name)
    padmet_file = all_study_padmet[study_name]
    if not os.path.isfile(sbml_file) and padmet_file:
        if verbose:
            print("Creating sbml from padmet for %s" %study_name)
        cmds = ["python3", padmet_utils_path + "/padmet_utils/connection/sbmlGenerator.py", "--padmet", padmet_file,
                "--output", sbml_file, "--sbml_lvl", "3", verbose]
        subprocess.call(cmds)


def convert_sbml_db(all_sbml_from_ortho):
    for sbml_file in all_sbml_from_ortho:
        if os.path.isfile(sbml_file):
            dict_file = "{0}_dict.csv".format(os.path.splitext(sbml_file)[0])
            if not os.path.exists(dict_file):
                cmds = ["python3", padmet_utils_path + "/padmet_utils/exploration/convert_sbml_db.py", "--mnx_rxn", mnx_rxn_path, "--sbml", sbml_file]
                db_ref = [line.split(":")[1] for line in subprocess.check_output(cmds, universal_newlines=True).splitlines() if line.startswith("Database")][0]
                if verbose:
                    print("%s: %s" %(os.path.basename(sbml_file), db_ref))
                if db_ref.lower() != "metacyc":
                    if verbose:
                        print("Creating id mapping file: %s" %dict_file)
                    cmds = ["python3",  padmet_utils_path+ "/padmet_utils/exploration/convert_sbml_db.py", "--mnx_rxn", mnx_rxn_path, "--mnx_cpd", mnx_cpd_path,
                            "--sbml", sbml_file, "--output", dict_file, "--db_out", "metacyc", verbose]
                    subprocess.call(cmds)


def create_draft(study_name):
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
        if os.path.exists(all_study_padmet[study_name]):
            if verbose:
                print("\tStarting from %s" %os.path.basename(all_study_padmet[study_name]))
            padmet_path = all_study_padmet[study_name]
            if os.path.exists(ortho_sbml_folder):
                cmd = "python3 {0}/padmet_utils/connection/sbml_to_padmet.py --padmetRef={1} --sbml={2} {3} --padmetSpec={4} --output={5} --source_tool={6} --source_category={7}".format(\
                padmet_utils_path, database_path, ortho_sbml_folder, verbose, padmet_path, output, source_tool, source_category)
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
                    "--padmetSpec", output, "--source_tool", source_tool, "--source_category", source_category, verbose]
        if os.path.exists(ortho_sbml_folder) and next(os.walk(ortho_sbml_folder))[2]:
            subprocess.call(cmds)
        else:
            if verbose:
                print("\t%s's folder is empty" %study_name)
            return


def modify_working_folder(working_folder):
    """
    Read this script.
    Search for the first line containing variable all_run_folder.
    Modify it to add the new working folder.
    Then rewrite the script.
    """
    with open(__file__, 'r') as aucome_script:
        aucome_lines = aucome_script.read().split('\n')

        for index, aucome_line in enumerate(aucome_lines):
            if '    all_run_folder = ' in aucome_line:
                aucome_line = '    all_run_folder = "{0}"'.format(working_folder)
                aucome_lines[index] = aucome_line
                break
        new_aucome_script_string = '\n'.join(aucome_lines)

    with open(__file__, 'w') as new_aucome_script:
        new_aucome_script.write(new_aucome_script_string)

    return


def installing_pwt(pwt_path, input_ptools_local_path):
    """
    Install silently Pathway-Tools in /programs.
    After running this function you need to source the bashrc.
    """
    ptools_local_path = '/root'
    if input_ptools_local_path:
        if os.path.isdir(input_ptools_local_path):
            ptools_local_path = input_ptools_local_path
        else:
            print(input_ptools_local_path + ' path does not exist, --ptools must be an existing path.')
            return

    cmd_chmods = ["chmod", "u+x", pwt_path]
    cmd_installs = [pwt_path, "--InstallDir", "/programs/pathway-tools", "--PTOOLS_LOCAL_PATH", ptools_local_path,
                    "--InstallDesktopShortcuts", "0", "--mode unattended"]
    cmd_echos = ["echo", 'export PATH="$PATH:/programs/pathway-tools:"', ">>", "~/.bashrc"]
    cmds = [cmd_chmods, cmd_installs, cmd_echos]
    for cmd in cmds:
        print(cmd)
        subprocess.call(cmd)
    print("Now you need to source your bash, run:")
    print("source ~/.bashrc")
    return


def uninstalling_pwt():
    """
    Uninstall Pathway-Tools and can delete ptools-local folder.
    """
    def ask_delete_ptools(ptools_path):
        yes_or_no = input('Delete ptools-local folder (y/n)?')
        if yes_or_no == 'y':
            subprocess.call(["rm",  "-r", ptools_path])
            print('Uninstallation of Pahtway-Tools and ptools-local done!')
            return
        elif yes_or_no == 'n':
            print('Uninstallation of Pathway-Tools done!.')
            return
        else:
            print('Wrong command')
            ask_delete_ptools(ptools_path)

    cmd_uninstall = ["/programs/pathway-tools/uninstall", "--mode unattended"]
    cmd_clean_bash = ["grep", "-v", """'export PATH="$PATH:/programs/pathway-tools:"'""", "~/.bashrc", ">", "~/temp.bashrc;", "mv", "~/temp.bashrc", "~/.bashrc"]
    cmds = [cmd_uninstall, cmd_clean_bash]

    ptools_path = mpwt.find_ptools_path()

    if os.path.isdir('/root/AIC-prefs'):
        cmd_delete_AIC_pref = ["rm", "-r" "/root/AIC-prefs"]
        cmds.append(cmd_delete_AIC_pref)

    for cmd in cmds:
        print(cmd)
        subprocess.call(cmd)

    ask_delete_ptools(ptools_path)
    return


def checking_genbank(genbank_file_name):
    """
    Check if there is a special character in the genbank file.
    If yes exit and print an error.

    Check gene ID in genbank to find too long gene ID or invalid character in gene ID.
    """
    invalid_characters = ['-', '|', '/', '(', ')', '\'', '=', '#', '*',
                '.', ':', '!', '+', '[', ']', ',', " "]
    if any(char in invalid_characters for char in genbank_file_name):
        print('Error in genbank file name: ' + genbank_file_name)
        print('Rename the file without:',genbank_file_name)

    # Path to the genbank file.
    genbank_path = studied_organisms_path + '/' + genbank_file_name + '/' + genbank_file_name + '.gbk'

    fix_dot_protein_seq = None
    invalid_gene_ids = []
    too_long_ids = []
    for record in SeqIO.parse(genbank_path, 'genbank'):
        for feature in record.features:
            if 'locus_tag' in feature.qualifiers:
                locus_tag = feature.qualifiers['locus_tag'][0]
                if any(char in invalid_characters for char in locus_tag):
                    if verbose:
                        invalid_gene_ids.append(feature)
                if len(locus_tag) >= 40:
                        too_long_ids.append(feature)
            if 'translation' in feature.qualifiers:
                if '.' in feature.qualifiers['translation'][0]:
                    fix_dot_protein_seq = True

    if len(invalid_gene_ids) > 0:
        print('Error of gene id in genbank ' + genbank_file_name + ', ' + str(len(invalid_gene_ids)) + ' genes have an invalid characters present: ' + ' '.join(invalid_characters) + '.')
    if len(too_long_ids) > 0:
        print('Error of gene id in genbank ' + genbank_file_name + ', ' + str(len(too_long_ids)) + ' genes have a gene id too long.')

    if len(invalid_gene_ids) > 0 or len(too_long_ids) > 0:
        print('Gene ID in ' + genbank_file_name + ' must be renamed.')
        fix_name = True
    else:
        fix_name = False

    if fix_dot_protein_seq:
        print('Dot in a protein sequence, Orthofinder will not work for this sequence. Dot will be deleted.')

    if fix_name or fix_dot_protein_seq:
        fix_genbank_file(genbank_file_name, fix_name, fix_dot_protein_seq)


def adapt_gene_id(gene_id, longest_gene_number_length):
    """
    Input: a gene ID like g_1 and the longest_gene_number_length (5 if you have a gene like g_12356).
    Return a new gene ID with a completion of 0 like: g_00001.
    """
    gene_prefix = '_'.join(gene_id.split('_')[:-1])
    gene_number = gene_id.split('_')[-1]

    new_gene_id = gene_prefix + '_' + gene_number.zfill(longest_gene_number_length)

    return new_gene_id


def fix_genbank_file(genbank_file_name, fix_name, fix_dot_protein_seq):
    # Path to the genbank file.
    genbank_path = studied_organisms_path + '/' + genbank_file_name + '/' + genbank_file_name + '.gbk'
    genbank_path_renamed = studied_organisms_path + '/' + genbank_file_name + '/' + genbank_file_name + '_original.gbk'

    if os.path.exists(genbank_path_renamed):
        print(genbank_file_name + ': Renaming has already been made on the data.')
        return

    # Create records that will be modified according to the issue: location or gene id.
    new_records = [record for record in SeqIO.parse(genbank_path, 'genbank')]

    # Use either the genbank accession or the genus + species as a prefix for new gene ID.
    try:
        new_prefix = new_records[0].annotations['accessions'][0] + '_' + str(new_records[0].annotations['sequence_version'])
    except:
        new_prefix = new_records[0].annotations['organism'].split(' ')[0][0] + '_' + new_records[0].annotations['organism'].split(' ')[1]

    if fix_name or fix_dot_protein_seq:
        # Dictionary wtih gene id as key and renamed id as value.
        feature_id_mappings = {}

        number_genes_genbanks = len([feature for record in SeqIO.parse(genbank_path, 'genbank') for feature in record.features if feature.type == 'gene'])
        gene_number = 1
        # Renamed ID: genbank file name + '_' + gene_position_number.
        # Max ID len is 39 for Pathway-Tools.
        for record in new_records:
            for feature in record.features:
                if 'locus_tag' in feature.qualifiers:
                    if fix_name:
                        feature_id = feature.qualifiers['locus_tag'][0]
                        if feature_id not in feature_id_mappings:
                            new_gene_id = new_prefix + '_' + str(gene_number)
                            new_feature_id = adapt_gene_id(new_gene_id, len(str(number_genes_genbanks)))
                            feature_id_mappings[feature_id] = new_feature_id
                            feature.qualifiers['locus_tag'][0] = new_feature_id
                            feature.qualifiers['old_locus_tag'] = feature_id
                            gene_number += 1
                        else:
                            feature.qualifiers['locus_tag'][0] = feature_id_mappings[feature_id]
                            feature.qualifiers['old_locus_tag'] = feature_id
                    if fix_dot_protein_seq:
                        if 'translation' in feature.qualifiers:
                            feature.qualifiers['translation'] = feature.qualifiers['translation'][0].replace('.', '')

        if fix_name:
            # Create a TSV mapping file with original and renamed ids.
            mapping_dic_path = studied_organisms_path + '/' + genbank_file_name + '/' + genbank_file_name + '_dict.csv'
            with open(mapping_dic_path, 'w') as csv_file:
                writer = csv.writer(csv_file, delimiter='\t')
                writer.writerow(["original_gene_id", "renamed_gene_id"])
                for key, value in list(feature_id_mappings.items()):
                    writer.writerow([key, value])

    # Create genbank with renamed id.
    new_genbank_path = studied_organisms_path + '/' + genbank_file_name + '/' + genbank_file_name + '_tmp.gbk'
    SeqIO.write(new_records, new_genbank_path, 'genbank')

    # Save original genbank.
    os.rename(genbank_path, genbank_path_renamed)

    # Rename renamed genbank to genbank used by the script.
    os.rename(new_genbank_path, genbank_path)

    if verbose:
        print(genbank_file_name + ' ids have been renamed.')


def orthogroup_to_sbml(dict_data):
    """
    dict_orthogroup: global var
    """
    #dict_data = {'study_name':'', 'o_compare_name': '', sbml_template':'', 'output':''}
    study_name = dict_data['study_name']
    to_compare_name = dict_data['to_compare_name']
    sbml_template = dict_data['sbml_template']
    output = dict_data['output']
    verbose = dict_data.get('verbose')
    dict_orthogroup = dict_data.get('ortho')
    if os.path.isfile(output):
        if verbose:
            print("*{0} is already created, skip".format(os.path.basename(output)))
        return
    if verbose:
        print("*Extracting orthology data to create sbml of {0} from {1}".format(study_name, to_compare_name))

    #k = gene_id from to_compare, v = list of genes id of study
    sub_dict_orth = {}
    for k in list(dict_orthogroup.values()):
        try:
            all_to_compare_genes = k[to_compare_name]
            all_study_genes = k[study_name]
            for to_compare_gene in all_to_compare_genes:
                try:
                    sub_dict_orth[to_compare_gene].update(all_study_genes)
                except KeyError:
                    sub_dict_orth[to_compare_gene] = set(all_study_genes)
        except KeyError:
            pass

    if not sub_dict_orth:
        if verbose:
            print("\t{0} and {1} don't share any ortholgue".format(study_name, to_compare_name))
        return
    print(sbml_template)
    reader = libsbml.SBMLReader()
    document_to_compare = reader.readSBML(sbml_template)
    for i in range(document_to_compare.getNumErrors()):
        print(document_to_compare.getError(i).getMessage())
    model_to_compare = document_to_compare.getModel()
    listOfReactions_with_genes = [rxn for rxn in model_to_compare.getListOfReactions()
                                  if sp.parseNotes(rxn).get("GENE_ASSOCIATION",[None])[0]]
    if verbose:
        print("\tSbml of {0} contains {1}/{2} reactions with genes assocation".format(to_compare_name, len(listOfReactions_with_genes), len(model_to_compare.getListOfReactions())))
    dict_rxn_ga = {}
    for rxn in listOfReactions_with_genes:
        ga = sp.parseNotes(rxn)['GENE_ASSOCIATION'][0]
        ga_for_gbr = re.sub(r" or " , "|", ga)
        ga_for_gbr = re.sub(r" and " , "&", ga_for_gbr)
        ga_for_gbr = re.sub(r"\s" , "", ga_for_gbr)
        #ga_for_gbr = "\"" + ga_for_gbr + "\""
        #ga_subsets = eval(subprocess.check_output("python3 grammar-boolean-rapsody.py "+ga_for_gbr, shell=True))
        if re.findall("\||&", ga_for_gbr):
            to_compare_ga_subsets = [_ for _ in gbr.compile_input(ga_for_gbr)]
        else:
            to_compare_ga_subsets = [(re.sub(r'\(|\)|\"', "", ga_for_gbr),)]     
        study_ga_subsets = []
        """
        to_compare_ga_subsets = [('a','c','d'),('c',)]
        sub_dict_orth = {'a':['a_a'],'c':['c_c'], 'd':['d_d']}
        """
        for to_compare_subset in to_compare_ga_subsets:
            study_subset = set()
            for gene in to_compare_subset:
                if gene in list(sub_dict_orth.keys()):
                    study_subset.update(sub_dict_orth[gene])
                else:
                    study_subset = set()
                    break
            if study_subset:
                """
                if verbose:
                    print("\t\t{0} == {1}".format(tuple(to_compare_subset), tuple(study_subset)))
                """
                study_ga_subsets.append(study_subset)
        if study_ga_subsets:
            study_ga = " or ".join(["("+" and ".join(subset)+")" for subset in study_ga_subsets])
            if verbose:
                print("\t\tAdding %s" %rxn.id)
                print("\t\tGENE_ASSOCIATION: %s" %(study_ga))
            dict_rxn_ga[rxn.id] = study_ga
    if not dict_rxn_ga:
        if verbose:
            print("\tNo reaction added from {0} to {1} because of missing orthologues".format(to_compare_name, study_name))
        return
    rxn_id_to_remove = set([rxn.id for rxn in model_to_compare.getListOfReactions()]).difference(list(dict_rxn_ga.keys()))
    if verbose:
        print("\tRemoving %s unused reactions" %len(rxn_id_to_remove))
    [model_to_compare.removeReaction(rxn_id) for rxn_id in rxn_id_to_remove]
    cpd_id_to_preserve = set()
    for rxn_id, study_ga in list(dict_rxn_ga.items()):
        rxn = model_to_compare.getElementBySId(rxn_id)
        #update notes
        notes_in_dict = sp.parseNotes(rxn)
        notes_in_dict["GENE_ASSOCIATION"] = [study_ga]
        notes = "<body xmlns=\"http://www.w3.org/1999/xhtml\">"
        for k,v_list in list(notes_in_dict.items()):
            for v in v_list:
                notes += "<p>"+k+": "+v+"</p>"
        notes += "</body>"
        rxn.setNotes(notes)
        cpd_in_rxn = set([p.getSpecies() for p in rxn.getListOfProducts()]).union(\
                         set([r.getSpecies() for r in rxn.getListOfReactants()]))
        cpd_id_to_preserve.update(cpd_in_rxn)
    all_species = [cpd.id for cpd in model_to_compare.getListOfSpecies()]
    [model_to_compare.removeSpecies(cpd_id ) for cpd_id in all_species if cpd_id not in cpd_id_to_preserve]
    new_id = os.path.basename(os.path.splitext(output)[0])    
    model_to_compare.setId(new_id)
    libsbml.writeSBMLToFile(document_to_compare, output)   
        

def create_config_file(config_file_path, run_id):
    config = configparser.RawConfigParser()
    config.add_section('DATABASE_PATHS')
    config.set('DATABASE_PATHS', 'database_ref_path', '/home/database/BIOCYC/METACYC/22.0_enhanced/metacyc_22.0_enhanced.padmet')
    config.set('DATABASE_PATHS', 'mnx_rxn_path', '/home/database/MNX/reac_xref.tsv')
    config.set('DATABASE_PATHS', 'mnx_cpd_path', '/home/database/MNX/chem_xref.tsv')
    config.add_section('PATHS_IN_RUN')
    config.set('PATHS_IN_RUN', 'run_id', run_id)
    config.set('PATHS_IN_RUN', 'studied_organisms_path', '%(run_id)s/studied_organisms')
    config.set('PATHS_IN_RUN', 'model_organisms_path', '%(run_id)s/model_organisms')
    config.set('PATHS_IN_RUN', 'orthology_based_path', '%(run_id)s/orthology_based')
    config.set('PATHS_IN_RUN', 'orthofinder_wd_path', '%(run_id)s/orthology_based/Orthofinder_WD')
    config.set('PATHS_IN_RUN', 'annotation_based_path', '%(run_id)s/annotation_based')
    config.set('PATHS_IN_RUN', 'pgdb_from_annotation_path', '%(annotation_based_path)s/PGDBs')
    config.set('PATHS_IN_RUN', 'padmet_from_annotation_path', '%(annotation_based_path)s/PADMETs')
    config.set('PATHS_IN_RUN', 'sbml_from_annotation_path', '%(annotation_based_path)s/SBMLs')
    config.set('PATHS_IN_RUN', 'networks_path', '%(run_id)s/networks')
    config.add_section('TOOL_PATHS')
    config.set('TOOL_PATHS', 'orthofinder_bin_path', '/programs/OrthoFinder-2.2.7/orthofinder')
    config.set('TOOL_PATHS', 'padmet_utils_path', '/programs/padmet-utils')
    config.add_section('VAR')
    config.set('VAR', 'study_from_annot_prefix', 'output_pathwaytools_')

    # Writing our configuration file to 'example.cfg'
    with open(config_file_path, 'w') as configfile:
        config.write(configfile)


def create_run(run_id):
    """
    create a run folder
    #todo, mkdir all tree
    """
    if os.path.isdir("{0}/{1}".format(all_run_folder, run_id)):
        print("Run '%s' already exist, remove this folder manually before" %run_id)
    else:
        print("creating Run %s" %run_id)
        os.mkdir("{0}/{1}".format(all_run_folder, run_id))
        all_folders = ["studied_organisms", "model_organisms", "networks", "orthology_based",\
                       "orthology_based/Orthofinder_WD", "annotation_based",\
                       "annotation_based/PGDBs", "annotation_based/PADMETs",\
                       "annotation_based/SBMLs", "analysis"]
        for folder in all_folders:
            print("creating folder {0}/{1}/{2}".format(all_run_folder, run_id, folder))
            os.mkdir("{0}/{1}/{2}".format(all_run_folder, run_id, folder))
        config_file_path = "{0}/{1}/config.txt".format(all_run_folder, run_id)
        create_config_file(config_file_path, run_id)


def get_version():
    '''
    Get version from Gitlab.
    Check internet connection using requests and eventlet timeout.
    '''
    reg_version = r'^\#+VERSION:([0-9.]*)#+'
    with eventlet.Timeout(2):
        try:
            response = requests.get(release_on_gitlab)
            first_line = response.text.split('\n')[0]
            version = re.match(reg_version,first_line).group(1)
        except eventlet.timeout.Timeout:
            version = None

    return version

if __name__ == "__main__":
    main()