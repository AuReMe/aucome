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
from padmet.utils import sbmlPlugin as sp

def main():
    args = docopt.docopt(__doc__)
    """
    args = {"--run_id":"test", "-v":True}
    all_run_folder = "/home/maite/Forge/docker/comparison_workspace"
    database_path = "/home/maite/Forge/docker/comparison_workspace/template/database/metacyc_22.0_enhanced.padmet"
    """
    run_id = args["--run_id"]
    global all_run_folder, database_path, studied_organisms_path, model_data_path, orthology_based_path, annotation_based_path,\
    sbml_from_annotation_path, networks_path, orthofinder_bin_path,\
    sbml_study_prefix, all_study_name, all_mode_name, dict_orthogroup, dir_path_gbr
    
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
    orthofinder_wd_path = "{0}/{1}".format(all_run_folder, config.get('RUN_PATHS','orthofinder_wd_path'))    
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
    dir_path_gbr = config.get('PADMET', 'dir_path_gbr')

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

    #for each faa, check if already in ortho_based
    if args["-o"]:
        #check if Orthofinder already run, if yes, get the last workdir

        for name, faa_path in all_study_faa.items():
            if not os.path.isfile("{0}/{1}.faa".format(orthofinder_wd_path, name)):
                if args["-v"]:
                    print("Copying {0}'s faa to {1}".format(name, orthofinder_wd_path))
                cmd = "cp {0} {1}/".format(faa_path, orthofinder_wd_path)
                subprocess.call(cmd, shell=True)
        for name, faa_path in all_model_faa.items():
            if not os.path.isfile("{0}/{1}.faa".format(orthofinder_wd_path, name)):
                if args["-v"]:
                    print("Copying {0}'s faa to {1}".format(name, orthofinder_wd_path))
                cmd = "cp {0} {1}/".format(faa_path, orthofinder_wd_path)
                subprocess.call(cmd, shell=True)
        if args["-v"]:
            print("Running Orthofinder")
        cmd = "{0} -f {1}".format(orthofinder_bin_path, orthofinder_wd_path)
        #check if folder "Result_..." in ortho_b_path/
        #2 diff command to run Orthofinder
        
        last_orthogroup_file = "/home/maite/Forge/docker/comparison_workspace/test_2/orthology_based/Orthofinder_WD/orthogroups.csv"
        if args["-v"]:
            print("Parsing Orthofinder output %s" %last_orthogroup_file)
        #k=orthogroup_id, v = {k = name, v = set of genes}
        dict_orthogroup = {}
        with open(last_orthogroup_file, 'r') as csvfile:
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
                os.makedirs(ortho_sbml_folder)
            all_to_compare = all_study_name.union(all_model_name) - set([study_name])
            all_dict_data = []
            for to_compare_name in all_to_compare:
                output = "{0}/{1}/output_orthofinder_{1}_from_{2}.sbml".format(orthology_based_path, study_name, to_compare_name)
                try:
                    sbml_template = all_model_sbml[to_compare_name]
                except KeyError:
                    sbml_template = all_study_sbml[to_compare_name]
                if sbml_template:
                    dict_data = {'study_name': study_name, 'to_compare_name': to_compare_name, 'sbml_template': sbml_template, 'output': output, 'verbose':True}
                    all_dict_data.append(dict_data)
            for dict_data in all_dict_data:
                orthogroup_to_sbml(dict_data)
            """
            p = Pool(processes=cpu_count())
            p.map_async(orthogroup_to_sbml, all_dict_data)
            p.close()
            p.join()
            """
                            
            

def orthogroup_to_sbml(dict_data):
    """
    dict_orthogroup: global var
    """
    #dict_data = {'study_name':'', 'o_compare_name': '', sbml_template':'', 'output':''}
    study_name = dict_data['study_name']
    to_compare_name = dict_data['to_compare_name']
    sbml_template = dict_data['sbml_template']
    output = dict_data['output']
    verbose = dict_data.get('verbose',False)
    #k = gene_id from to_compare, v = list of genes id of study
    sub_dict_orth = {}
    for k in dict_orthogroup.values():
        all_to_compare_genes = k[to_compare_name]
        all_study_genes = k[study_name]
        if all_to_compare_genes and all_study_genes:
            for to_compare_gene in all_to_compare_genes:
                try:
                    sub_dict_orth[to_compare_gene].update(all_study_genes)
                except KeyError:
                    sub_dict_orth[to_compare_gene] = set(all_study_genes)
                    
    if not sub_dict_orth:
        if verbose:
            print("{0} and {1} don't share any ortholgue".format(study_name, to_compare_name))
        return

    reader = libsbml.SBMLReader()
    document_to_compare = reader.readSBML(sbml_template)
    for i in range(document_to_compare.getNumErrors()):
        print (document_to_compare.getError(i).getMessage())
    model_to_compare = document_to_compare.getModel()
    listOfReactions_with_genes = [rxn for rxn in model_to_compare.getListOfReactions()
                                  if sp.parseNotes(rxn).get("GENE_ASSOCIATION",[None])[0]]
    if verbose:
        print("SBML Model contains %s/%s reactions with genes assocation" %(len(listOfReactions_with_genes), len(model_to_compare.getListOfReactions())))
    dict_rxn_ga = {}
    for rxn in listOfReactions_with_genes:
        ga = sp.parseNotes(rxn)['GENE_ASSOCIATION'][0]
        ga_for_gbr = re.sub(r" or " , "|", ga)
        ga_for_gbr = re.sub(r" and " , "&", ga_for_gbr)
        ga_for_gbr = re.sub(r"\s" , "", ga_for_gbr)
        ga_for_gbr = "\"" + ga_for_gbr + "\""
        #ga_subsets = eval(subprocess.check_output("python3 grammar-boolean-rapsody.py "+ga_for_gbr, shell=True))
        if re.findall("\||&", ga_for_gbr):
            to_compare_ga_subsets = eval(subprocess.check_output("python3 "+dir_path_gbr+" "+ga_for_gbr, shell=True))
        else:
            to_compare_ga_subsets = [ga_for_gbr]            
        study_ga_subsets = []
        """
        to_compare_ga_subsets = [('a','c','d'),('c',)]
        sub_dict_orth = {'a':['a_a'],'c':['c_c'], 'd':['d_d']}
        """
        for to_compare_subset in to_compare_ga_subsets:
            study_subset = set()
            for gene in to_compare_subset:
                if gene in sub_dict_orth.keys():
                    study_subset.update(sub_dict_orth[gene])
                else:
                    study_subset = set()
                    break
            if study_subset:
                if verbose:
                    print("{0} == {1}".format(tuple(to_compare_subset), tuple(study_subset)))
                study_ga_subsets.append(study_subset)
        if study_ga_subsets:
            study_ga = " or ".join(["("+" and ".join(subset)+")" for subset in study_ga_subsets])
            if verbose:
                print("Adding %s" %rxn.id)
                print("GENE_ASSOCIATION: %s" %(study_ga))
            dict_rxn_ga[rxn.id] = study_ga
    if not dict_rxn_ga:
        if verbose:
            print("No reaction added because of missing orthologues")
        return
    rxn_id_to_remove = set([rxn.id for rxn in model_to_compare.getListOfReactions()]).difference(dict_rxn_ga.keys())
    if verbose:
        print("Removing unused reactions" %len(rxn_id_to_remove))
    [model_to_compare.removeReaction(rxn_id) for rxn_id in rxn_id_to_remove]
    cpd_id_to_preseve = set()
    for rxn_id in dict_rxn_ga.keys():
        rxn = model_to_compare.getElementBySId(rxn_id)
        cpd_in_rxn = set([p.getSpecies() for p in rxn.getListOfProducts()]).union(\
                         set([r.getSpecies() for r in rxn.getListOfReactants()]))
        cpd_id_to_preseve.update(cpd_in_rxn)
    [model_to_compare.removeSpecies(cpd.id ) for cpd in model_to_compare.getListOfSpecies() if cpd.id not in cpd_id_to_preseve]
        
    for rxn_id, gene_assoc in  dict_rxn_ga.items():
        
                

    
def create_config_file(config_file_path, run_id):
    config = configparser.RawConfigParser()
    config.add_section('RUN_PATHS')
    config.set('RUN_PATHS', 'run_id', run_id)
    config.set('RUN_PATHS', 'database_path', '/home/data/database/BIOCYC/Metacyc/22.0/metacyc_22.0_enhanced.padmet')
    config.set('RUN_PATHS', 'studied_organisms_path', '%(run_id)s/studied_organisms')
    config.set('RUN_PATHS', 'model_organisms_path', '%(run_id)s/model_organisms')
    config.set('RUN_PATHS', 'orthology_based_path', '%(run_id)s/orthology_based')
    config.set('RUN_PATHS', 'orthofinder_wd_path', '%(run_id)s/orthology_based/Orthofinder_WD')
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
    config.set('PADMET', 'dir_path_gbr', '%(run_id)s/connection/grammar-boolean-rapsody.py')
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