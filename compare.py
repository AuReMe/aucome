# -*- coding: utf-8 -*-
"""
Created on Thu Apr 26 16:32:23 2018

@author: maite

1./ get all fasta in studied_orgnisms/ in dict_faa_paths['study']
ex dict_faa_paths['study'] = [studied_organisms_path/study_1/faa_study_name, ...]
get all fasta in model_data/
ex dict_faa_paths['model'] = [model_data_path/model_1/faa_model_name, ...]

usage:
    compare.py --run_id=DIR [-i] [-o] [-d] [-v] [--log=FILE]

options:
    -h --help     Show help.
    --run_id=FILE    pathname to the comparison workspace
    -o    Run Orthofinder
    -p    Run pantograph
    -d    Run inparanoid, omlc, pantograph and merge all networks
    -v    Verbose

"""
import docopt
import os
import subprocess
import itertools
from multiprocessing import Pool, cpu_count
import time
import configparser
import re
import libsbml
import csv

def main():
    args = docopt.docopt(__doc__)
    """
    args = {"--run_id":"test_2", "-v":True}
    all_run_folder = "/home/maite/Forge/docker/comparison_workspace"
    database_path = "/home/maite/Forge/docker/comparison_workspace/template/database/metacyc_22.0_enhanced.padmet"
    """
    run_id = args["--run_id"]
    global all_run_folder, database_path, studied_organisms_path, model_data_path, orthology_based_path, annotation_based_path,\
    sbml_from_annotation_path, networks_path, orthofinder_bin_path,\
    sbml_study_prefix, all_study_name, all_mode_name
    
    all_run_folder = "/shared"
    config_path = "{0}/{1}/config.txt".format(all_run_folder, run_id)
    config = configparser.ConfigParser()
    config.read(config_path)
    #config.read('/home/maite/Forge/docker/comparison_workspace/template/config.txt')
    #RUN_PATHS
    database_path = config.get('RUN_PATHS',"database_path")
    studied_organisms_path = "{0}/{1}".format(all_run_folder, config.get('RUN_PATHS','studied_organisms_path'))
    model_organisms_path = "{0}/{1}".format(all_run_folder, config.get('RUN_PATHS','model_organisms_path'))
    orthology_based_path = "{0}/{1}".format(all_run_folder, config.get('RUN_PATHS','orthology_based_path'))
    annotation_based_path = "{0}/{1}".format(all_run_folder, config.get('RUN_PATHS','annotation_based_path'))
    padmet_from_annotation_path = "{0}/{1}".format(all_run_folder, config.get('RUN_PATHS','padmet_from_annotation_path'))
    sbml_from_annotation_path = "{0}/{1}".format(all_run_folder, config.get('RUN_PATHS','sbml_from_annotation_path'))
    networks_path = "{0}/{1}".format(all_run_folder, config.get('RUN_PATHS','networks_path'))
    #ORTHO
    orthofinder_bin_path = config.get('ORTHOLOGY','orthofinder_bin_path')
    #ANNOT
    sbml_study_from_annot_prefix = config.get('ANNOTATION','sbml_study_from_annot_prefix')
    #PADMET
    padmet_utils_path = config.get('PADMET', 'padmet_utils_path')


    #create dict for ortho data
    all_study_name = set(os.walk(studied_organisms_path).next()[1])
    all_model_name = set(os.walk(model_organisms_path).next()[1])
    
    #check for each study if exist PGDB folder in PGDBs folder, if missing RUN ptools    


    all_study_gbk = dict([(study_name, "{0}/{1}/{1}.gbk".format(studied_organisms_path, study_name))
                          if os.path.isfile("{0}/{1}/{1}.gbk".format(studied_organisms_path, study_name))
                          else (study_name, None)
                          for study_name in all_study_name])

    #create Faa from gbk if no faa found    
    for study_name in all_study_name:
        faa_path = "{0}/{1}/{1}.faa".format(studied_organisms_path, study_name)
        gbk_file = all_study_gbk[study_name]
        if not os.path.isfile(faa_path) and gbk_file:
            if args["-v"]:
                print("Creating faa from gbk for %s" %study_name)
            cmd = "python {0}/connection/gbk_to_faa.py --gbk={1} --output={2}".format(padmet_utils_path, gbk_file, faa_path)
            subprocess.call(cmd, shell=True)

    #k = folder_name in studied_org_path, v = path to faa in this folder, faa name should be folder_name.faa
    all_study_faa = dict([(study_name, "{0}/{1}/{1}.faa".format(studied_organisms_path, study_name))
                          if os.path.isfile("{0}/{1}/{1}.faa".format(studied_organisms_path, study_name))
                          else (study_name, None)
                          for study_name in all_study_name])

    #sbml of study are obtained from annotation, they should be in sbml_from_annotation_path
    #k = study_name (== folder_name in studied_org_path or obtained from sbml name), v = path to sbml, sbml_study_prefi+study_name+.sbml
    all_study_sbml = dict([(study_name, "{0}/{1}{2}.sbml".format(sbml_from_annotation_path, sbml_study_from_annot_prefix, study_name))
                           if os.path.isfile("{0}/{1}{2}.sbml".format(sbml_from_annotation_path, sbml_study_from_annot_prefix, study_name))
                           else (study_name,None)
                           for study_name in all_study_name])

    #k = folder_name in model_organisms_path, v = path to faa in this folder, faa name should be folder_name.faa
    all_model_gbk = dict([(model_name, "{0}/{1}/{1}.gbk".format(model_organisms_path, model_name))
                          if os.path.isfile("{0}/{1}/{1}.gbk".format(model_organisms_path, model_name))
                          else (model_name, None)
                          for model_name in all_model_name])

    for model_name in all_model_name:
        faa_path = "{0}/{1}/{1}.faa".format(model_organisms_path, model_name)
        gbk_file = all_model_gbk[model_name]
        if not os.path.isfile(faa_path) and gbk_file:
            if args["-v"]:
                print("Creating faa from gbk for %s" %model_name)
            cmd = "python {0}/connection/gbk_to_faa.py --gbk={1} --output={2}".format(padmet_utils_path, gbk_file, faa_path)
            subprocess.call(cmd, shell=True)

    #k = folder_name in model_organisms_path, v = path to faa in this folder, faa name should be folder_name.faa
    all_model_faa = dict([(model_name, "{0}/{1}/{1}.faa".format(model_organisms_path, model_name))
                          if os.path.isfile("{0}/{1}/{1}.faa".format(model_organisms_path, model_name))
                          else (model_name, None)
                          for model_name in all_model_name])

    #k = folder_name in model_organisms_path, v = path to sbml in this folder, sbml name should be folder_name.sbml
    all_model_sbml = dict([(model_name, "{0}/{1}/{1}.sbml".format(model_organisms_path, model_name))
                          if os.path.isfile("{0}/{1}/{1}.sbml".format(model_organisms_path, model_name))
                          else (model_name, None)
                          for model_name in all_model_name])


    if args["-v"]:
        print("Input summary:")
        print("* %s Studied organims:" %(len(all_study_name)))
        for study_name in all_study_name:
            print("%s:" %study_name)
            if all_study_gbk[study_name]:
                print("\tGBK: OK")
            else:
                print("\t[WARNING] No GBK found, should be in {1}/{0}/{0}.gbk".format(study_name, studied_organisms_path))
            if all_study_faa[study_name]:
                print("\tFAA: OK")
            else:
                print("\t[WARNING] No FAA found, should be in {1}/{0}/{0}.faa".format(study_name, studied_organisms_path))
            if all_study_sbml[study_name]:
                print("\tSBML: OK")
            else:
                print("\t[WARNING] No SBML found, should be in {1}/{2}{0}.sbml".format(study_name, sbml_from_annotation_path, sbml_study_from_annot_prefix))
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
    exit()
    #for each faa, check if already in ortho_based
    if args["-o"]:
        #check if Orthofinder already run, if yes, get the last workdir

        for name, faa_path in all_study_faa.items():
            if not os.path.isfile("{0}/{1}.faa".format(orthology_based_path, name)):
                if args["-v"]:
                    print("Copying {0}'s faa to {1}".format(name, orthology_based_path))
                cmd = "cp {0} {1}/".format(faa_path, orthology_based_path)
                #subprocess.call(cmd, shell=True)
        for name, faa_path in all_model_faa.items():
            if not os.path.isfile("{0}/{1}.faa".format(orthology_based_path, name)):
                if args["-v"]:
                    print("Copying {0}'s faa to {1}".format(name, orthology_based_path))
                cmd = "cp {0} {1}/".format(faa_path, orthology_based_path)
                #subprocess.call(cmd, shell=True)
        if args["-v"]:
            print("Running Orthofinder")
        cmd = "{0} -f {1}".format(orthofinder_bin_path, orthology_based_path)
        #check if folder "Result_..." in ortho_b_path/
        #2 diff command to run Orthofinder
        
        last_orthogroup_file = "/home/maite/Forge/docker/comparison_workspace/OrthoFinder-2.2.6/ExampleData/Results_May23/Orthogroups.csv"
        if args["-v"]:
            print("Parsing Orthofinder output %s" %last_orthogroup_file)
        #k=orthogroup_id, v = {k = name, v = set of genes}
        dict_orthogroup = {}
        with open(orthofinder_output, 'r') as csvfile:
            reader = csv.DictReader(csvfile, delimiter = "\t")
            for row in reader:
                orth_id = row['']
                row.pop('')
                new_dict = dict([(name, set(genes.split(","))) for (name, genes) in row.items()])
                dict_orthogroup[orth_id] = new_dict

        if args["-v"]:
            print("Start multiprocess sbml creation...")
        for study_name in all_study_name:
            ortho_sbml_folder = "{0}/{1}".format(orthology_based_path, study_name)
            if not os.path.exists(ortho_sbml_folder):
                if args["-v"]:
                    print("Creating folder %s" %ortho_sbml_folder)
                os.makedirs(full_orth_run_folder)
            for to_compare in all_study_name.union(all_model_name) - set([study_name]):
                print("ok")
                
            

def parse_orthogroups(orthofinder_output):
    orthofinder_output = "/home/maite/Forge/docker/comparison_workspace/OrthoFinder-2.2.6/ExampleData/Results_May23/Orthogroups.csv"
    for study_name in all_study_name:
        sbml_template = "/home/maite/Forge/docker/comparison_workspace/template/model_organisms/Nannochloropsis_salina/Nannochloropsis_salina.sbml"
    with open(orthofinder_output, 'r') as csvfile:
        reader = csv.DictReader(csvfile, delimiter = "\t")
        print reader.fieldnames
        a = [row for row in reader]

    ga = "( Cre03.g144627.t1.1 and ( Cre10.g446100.t1.2 or Cre05.g243050.t1.2 or Cre01.g066552.t1.1 ) )"
    dir_path_gbr = "/home/maite/padmet-tools/padmet-utils/connection/grammar-boolean-rapsody.py"
    ga_for_gbr = re.sub(r" or " , "|", ga)
    ga_for_gbr = re.sub(r" and " , "&", ga_for_gbr)
    ga_for_gbr = re.sub(r"\s" , "", ga_for_gbr)
    ga_for_gbr = "\"" + ga_for_gbr + "\""
    #for edi test only
    #ga_subsets = eval(subprocess.check_output("python3 grammar-boolean-rapsody.py "+ga_for_gbr, shell=True))
    ga_subsets = eval(subprocess.check_output("python3 "+dir_path_gbr+" "+ga_for_gbr, shell=True))


    
def create_config_file(config_file_path, run_id):
    config = configparser.RawConfigParser()
    config.add_section('RUN_PATHS')
    config.set('RUN_PATHS', 'run_id', run_id)
    config.set('RUN_PATHS', 'database_path', '/home/data/database/BIOCYC/Metacyc/22.0/metacyc_22.0_enhanced.padmet')
    config.set('RUN_PATHS', 'studied_organisms_path', '%(run_id)s/studied_organisms')
    config.set('RUN_PATHS', 'model_organisms_path', '%(run_id)s/model_organisms')
    config.set('RUN_PATHS', 'orthology_based_path', '%(run_id)s/orthology_based')
    config.set('RUN_PATHS', 'annotation_based_path', '%(run_id)s/annotation_based')
    config.set('RUN_PATHS', 'padmet_from_annotation_path', '%(annotation_based_path)s/PADMETs')
    config.set('RUN_PATHS', 'sbml_from_annotation_path', '%(annotation_based_path)s/SBMLs')
    config.set('RUN_PATHS', 'networks_path', '%(run_id)s/networks')
    config.add_section('ORTHOLOGY')
    config.set('ORTHOLOGY', 'orthofinder_bin_path', '/programs/OrthoFinder-2.2.6/orthofinder')
    config.add_section('ANNOTATION')
    config.set('ANNOTATION', 'sbml_study_from_annot_prefix', 'output_pathwaytools_')
    config.add_section('PADMET')
    config.set('PADMET', 'padmet_utils_path', '/programs/padmet-utils')

    # Writing our configuration file to 'example.cfg'
    with open(config_file_path, 'wb') as configfile:
        config.write(configfile)

def create_run(run_id):
    """
    create a run folder
    #todo, mkdir all tree
    """
        



if __name__ == "__main__":
    main()