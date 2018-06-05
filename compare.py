#!/usr/bin/env python
# -*- coding: utf-8 -*-
#!/usr/bin/env python
"""
Created on Thu Apr 26 16:32:23 2018

@author: maite

1./ get all fasta in studied_orgnisms/ in dict_faa_paths['study']
ex dict_faa_paths['study'] = [studied_organisms_path/study_1/faa_study_name, ...]
get all fasta in model_data/
ex dict_faa_paths['model'] = [model_data_path/model_1/faa_model_name, ...]

usage:
    compare.py --init=ID [-v]
    compare.py --run=DIR [-c] [-o] [-p] [-d] [-v] [--log=FILE]
    compare.py -R

options:
    -h --help     Show help.
    -R     Open access from container.
    --run=ID    pathname to the comparison workspace
    -c    Check inputs validity
    -o    Run Orthofinder
    -p    Run Pathway-Tools
    -d    Run Orthofinder, Pathway and merge all networks
    -v     Verbose.

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
from padmet.classes import PadmetSpec, PadmetRef
#import mpwt

def main():
    args = docopt.docopt(__doc__)
    """
    args = {"--run":"algues", "-v":True}
    all_run_folder = "/home/maite/Forge/docker/comparison_workspace"
    database_path = "/home/maite/Forge/docker/comparison_workspace/template/database/metacyc_22.0_enhanced.padmet"
    padmet_utils_path = "/home/maite/padmet-tools/padmet-utils"
    """
    global all_run_folder, database_path, studied_organisms_path, model_data_path, orthology_based_path, annotation_based_path,\
    sbml_from_annotation_path, networks_path, orthofinder_bin_path,\
    sbml_study_prefix, all_study_name, all_mode_name, dict_orthogroup, padmet_utils_path

    all_run_folder = "/shared"
    #add permission to all folder in all_run_folder, usefull because all cmd exec from container are root based
    if args['-R']:
        cmd = "chmod -R 777 %s" %all_run_folder
        subprocess.call(cmd, shell=True)
        return

    if args["--run"]:
        run_id = args["--run"]
    else:
        run_id = args["--init"]
        
    config_file_path = "{0}/{1}/config.txt".format(all_run_folder, run_id)
    if args["--init"]:
        create_run(run_id)
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
    all_study_name = set(os.walk(studied_organisms_path).next()[1])
    all_model_name = set(os.walk(model_organisms_path).next()[1])

    if args["-p"]:
        #check for each study if exist PGDB folder in PGDBs folder, if missing RUN ptools
        chronoDepart = time.time()
        mpwt.multiprocess_pwt(input_folder=studied_organisms_path, output_folder=pgdb_from_annotation_path, dat_extraction=True)
        chrono = (time.time() - chronoDepart)
        partie_entiere, partie_decimale = str(chrono).split('.')
        chrono = ".".join([partie_entiere, partie_decimale[:3]])
        if verbose:
            print("Pathwaytools done in: %ss" %chrono)
    #PGDB, padmet, sbml
    all_study_pgdb = dict([(study_name, "{0}/{1}".format(pgdb_from_annotation_path, study_name))
                          if os.path.isdir("{0}/{1}".format(pgdb_from_annotation_path, study_name))
                          else (study_name, '')
                          for study_name in all_study_name])

    if args["-c"]:
        for study_name in all_study_name:
            padmet_file = "{0}/{1}{2}.padmet".format(padmet_from_annotation_path, study_from_annot_prefix, study_name)
            pgdb_folder = all_study_pgdb[study_name]
            if not os.path.isfile(padmet_file) and pgdb_folder:
                if verbose:
                    print("Creating padmet from pgdb for %s" %study_name)
                cmd = "python {0}/connection/pgdb_to_padmet.py --output={1} --directory={2} --padmetRef={3} --source=genome -g {4}".\
                format(padmet_utils_path, padmet_file, pgdb_folder, database_path, verbose)
                subprocess.call(cmd, shell=True)
    all_study_padmet = dict([(study_name, "{0}/{1}{2}.padmet".format(padmet_from_annotation_path, study_from_annot_prefix, study_name))
                          if os.path.isfile("{0}/{1}{2}.padmet".format(padmet_from_annotation_path, study_from_annot_prefix, study_name))
                          else (study_name, '')
                          for study_name in all_study_name])
    if args["-c"]:
        for study_name in all_study_name:
            sbml_file = "{0}/{1}{2}.sbml".format(sbml_from_annotation_path, study_from_annot_prefix, study_name)
            padmet_file = all_study_padmet[study_name]
            if not os.path.isfile(sbml_file) and padmet_file:
                if verbose:
                    print("Creating sbml from padmet for %s" %study_name)
                cmd = "python {0}/connection/sbmlGenerator.py --padmet={1} --output={2} --sbml_lvl=2 {3}".format(padmet_utils_path, padmet_file, sbml_file, verbose)
                subprocess.call(cmd, shell=True)

    #sbml of study are obtained from annotation, they should be in sbml_from_annotation_path
    #k = study_name (== folder_name in studied_org_path or obtained from sbml name), v = path to sbml, sbml_study_prefi+study_name+.sbml
    all_study_sbml = dict([(study_name, "{0}/{1}{2}.sbml".format(sbml_from_annotation_path, study_from_annot_prefix, study_name))
                           if os.path.isfile("{0}/{1}{2}.sbml".format(sbml_from_annotation_path, study_from_annot_prefix, study_name))
                           else (study_name, '')
                           for study_name in all_study_name])

    all_study_gbk = dict([(study_name, "{0}/{1}/{1}.gbk".format(studied_organisms_path, study_name))
                          if os.path.isfile("{0}/{1}/{1}.gbk".format(studied_organisms_path, study_name))
                          else (study_name, '')
                          for study_name in all_study_name])

    #create Faa from gbk if no faa found
    if args["-c"]:
        for study_name in all_study_name:
            faa_path = "{0}/{1}/{1}.faa".format(studied_organisms_path, study_name)
            gbk_file = all_study_gbk[study_name]
            if not os.path.isfile(faa_path) and gbk_file:
                if verbose:
                    print("Creating faa from gbk for %s" %study_name)
                cmd = "python {0}/connection/gbk_to_faa.py --gbk={1} --output={2}".format(padmet_utils_path, gbk_file, faa_path)
                subprocess.call(cmd, shell=True)

    #k = folder_name in studied_org_path, v = path to faa in this folder, faa name should be folder_name.faa
    all_study_faa = dict([(study_name, "{0}/{1}/{1}.faa".format(studied_organisms_path, study_name))
                          if os.path.isfile("{0}/{1}/{1}.faa".format(studied_organisms_path, study_name))
                          else (study_name, '')
                          for study_name in all_study_name])

    #k = folder_name in model_organisms_path, v = path to faa in this folder, faa name should be folder_name.faa
    all_model_gbk = dict([(model_name, "{0}/{1}/{1}.gbk".format(model_organisms_path, model_name))
                          if os.path.isfile("{0}/{1}/{1}.gbk".format(model_organisms_path, model_name))
                          else (model_name, '')
                          for model_name in all_model_name])

    if args["-c"]:
        for model_name in all_model_name:
            faa_path = "{0}/{1}/{1}.faa".format(model_organisms_path, model_name)
            gbk_file = all_model_gbk[model_name]
            if not os.path.isfile(faa_path) and gbk_file:
                if verbose:
                    print("Creating faa from gbk for %s" %model_name)
                cmd = "python {0}/connection/gbk_to_faa.py --gbk={1} --output={2}".format(padmet_utils_path, gbk_file, faa_path)
                subprocess.call(cmd, shell=True)

    #k = folder_name in model_organisms_path, v = path to faa in this folder, faa name should be folder_name.faa
    all_model_faa = dict([(model_name, "{0}/{1}/{1}.faa".format(model_organisms_path, model_name))
                          if os.path.isfile("{0}/{1}/{1}.faa".format(model_organisms_path, model_name))
                          else (model_name, '')
                          for model_name in all_model_name])

    #k = folder_name in model_organisms_path, v = path to sbml in this folder, sbml name should be folder_name.sbml
    all_model_sbml = dict([(model_name, "{0}/{1}/{1}.sbml".format(model_organisms_path, model_name))
                          if os.path.isfile("{0}/{1}/{1}.sbml".format(model_organisms_path, model_name))
                          else (model_name, '')
                          for model_name in all_model_name])


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
                print("\t[WARNING] No PADMET found, should be in {1}/{2}{0}.padmet".format(study_name, sbml_from_annotation_path, study_from_annot_prefix))
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

    #for each faa, check if already in ortho_based
    if args["-o"]:
        #check if Orthofinder already run, if yes, get the last workdir
        try:
            last_orthogroup_file = max(["%s/%s" %(x[0], 'Orthogroups.csv') for x in os.walk(orthofinder_wd_path) if 'Orthogroups.csv' in x[2]])
        except ValueError:
            if verbose:
                print("Enable to find file Orthogroups.csv in {1}, need to run Orthofinder...".format(orthofinder_wd_path))
            for name, faa_path in all_study_faa.items():
                if not os.path.isfile("{0}/{1}.faa".format(orthofinder_wd_path, name)):
                    if verbose:
                        print("Copying {0}'s faa to {1}".format(name, orthofinder_wd_path))
                    cmd = "cp {0} {1}/".format(faa_path, orthofinder_wd_path)
                    subprocess.call(cmd, shell=True)
            for name, faa_path in all_model_faa.items():
                if not os.path.isfile("{0}/{1}.faa".format(orthofinder_wd_path, name)):
                    if verbose:
                        print("Copying {0}'s faa to {1}".format(name, orthofinder_wd_path))
                    cmd = "cp {0} {1}/".format(faa_path, orthofinder_wd_path)
                    subprocess.call(cmd, shell=True)
            nb_cpu_to_use = cpu_count()-1
            if verbose:
                print("Running Orthofinder on %s cpu" %nb_cpu_to_use)

            chronoDepart = time.time()
            cmd = "{0} -f {1} -t {2}".format(orthofinder_bin_path, orthofinder_wd_path, nb_cpu_to_use)
            subprocess.call(cmd, shell=True)
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
                new_dict = dict([(name, set(genes.split(","))) for (name, genes) in row.items() if genes])
                dict_orthogroup[orth_id] = new_dict

        if verbose:
            print("Start multiprocess sbml creation...")
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
                    dict_data = {'study_name': study_name, 'to_compare_name': to_compare_name, 'sbml_template': sbml_template, 'output': output, 'verbose':True}
                    all_dict_data.append(dict_data)

            chronoDepart = time.time()
            p = Pool(processes=cpu_count())
            p.map_async(orthogroup_to_sbml, all_dict_data)
            p.close()
            p.join()
            """
            for dict_data in all_dict_data:
                orthogroup_to_sbml(dict_data)
            """
            chrono = (time.time() - chronoDepart)
            partie_entiere, partie_decimale = str(chrono).split('.')
            chrono = ".".join([partie_entiere, partie_decimale[:3]])
            if verbose:
                print("Orthofinder output parsed in: %ss" %chrono)
            #check database, mapping to metacyc ???
        all_sbml_from_ortho = [dict_data['output'] for dict_data in all_dict_data]
        for sbml_file in all_sbml_from_ortho:
            if os.path.isfile(sbml_file):
                cmd = "python {0}/exploration/convert_sbml_db.py --mnx_rxn={1} --sbml={2}".format(padmet_utils_path, mnx_rxn_path, sbml_file)
                db_ref = [line.split(":")[1] for line in subprocess.check_output(cmd, shell=True).splitlines() if line.startswith("Database")][0]
                if verbose:
                    print("%s: %s" %(sbml_file, db_ref))
                if db_ref.lower() != "metacyc":
                    dict_file = "{0}_dict.csv".format(os.path.splitext(sbml_file)[0])
                    if verbose:
                        print("Creating id mapping file: %s" %dict_file)
                    cmd = "python {0}/exploration/convert_sbml_db.py --mnx_rxn={1} --mnx_cpd={2} --sbml={3} --output={4} --db_out='metacyc' {5}".format(\
                    padmet_utils_path, mnx_rxn_path, mnx_cpd_path, sbml_file, dict_file, verbose)
                    
                    

    if args["-d"]:
        padmetRef = PadmetRef(database_path)
        for study_name in all_study_name:
            output = "{0}/{1}.padmet".format(networks_path, study_name)
            if os.path.exists(output):
                if verbose:
                    print("%s already exist, skip" %os.path.basename(output))
                    pass
            else:
                if verbose:
                    print("Creating %s" %os.path.basename(output))
                if os.path.exists(all_study_padmet[study_name]):
                    if verbose:
                        print("\tStarting from %s" %os.path.basename(all_study_padmet[study_name]))
                    padmet_path = all_study_padmet[study_name]
                    padmet = PadmetSpec(padmet_path)
                else:
                    if verbose:
                        print("\tStarting from an empty PADMET")
                    padmet = PadmetSpec()
                    padmet.setInfo(padmetRef)
                    padmet.setPolicy(padmetRef)                
                ortho_sbml_folder = "{0}/{1}".format(orthology_based_path, study_name)
                if os.path.exists(ortho_sbml_folder) and os.walk(ortho_sbml_folder).next()[2]:
                    cmd = ""
                    subprocess.call(cmd, shell=True)
                else:
                    if verbose:
                        print("\t%s's folder is empty" %study_name)
                    pass
            

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
    if verbose:
        print("%s:" %study_name)
    if os.path.isfile(output):
        if verbose:
            print("\t%s is already created, skip" %os.path.basename(output))
        return
    if verbose:
        print("\tExtracting orthology data to create sbml of {0} from {1}".format(study_name, to_compare_name))
    #k = gene_id from to_compare, v = list of genes id of study
    sub_dict_orth = {}
    for k in dict_orthogroup.values():
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

    reader = libsbml.SBMLReader()
    document_to_compare = reader.readSBML(sbml_template)
    for i in range(document_to_compare.getNumErrors()):
        print (document_to_compare.getError(i).getMessage())
    model_to_compare = document_to_compare.getModel()
    listOfReactions_with_genes = [rxn for rxn in model_to_compare.getListOfReactions()
                                  if sp.parseNotes(rxn).get("GENE_ASSOCIATION",[None])[0]]
    if verbose:
        print("\tSBML Model contains %s/%s reactions with genes assocation" %(len(listOfReactions_with_genes), len(model_to_compare.getListOfReactions())))
    dict_rxn_ga = {}
    for rxn in listOfReactions_with_genes:
        ga = sp.parseNotes(rxn)['GENE_ASSOCIATION'][0]
        ga_for_gbr = re.sub(r" or " , "|", ga)
        ga_for_gbr = re.sub(r" and " , "&", ga_for_gbr)
        ga_for_gbr = re.sub(r"\s" , "", ga_for_gbr)
        ga_for_gbr = "\"" + ga_for_gbr + "\""
        #ga_subsets = eval(subprocess.check_output("python3 grammar-boolean-rapsody.py "+ga_for_gbr, shell=True))
        if re.findall("\||&", ga_for_gbr):
            cmd = "python3 {0}/connection/grammar-boolean-rapsody.py {1}".format(padmet_utils_path, ga_for_gbr)
            to_compare_ga_subsets = eval(subprocess.check_output(cmd, shell=True))
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
    rxn_id_to_remove = set([rxn.id for rxn in model_to_compare.getListOfReactions()]).difference(dict_rxn_ga.keys())
    if verbose:
        print("\tRemoving %s unused reactions" %len(rxn_id_to_remove))
    [model_to_compare.removeReaction(rxn_id) for rxn_id in rxn_id_to_remove]
    cpd_id_to_preserve = set()
    for rxn_id, study_ga in dict_rxn_ga.items():
        rxn = model_to_compare.getElementBySId(rxn_id)
        #update notes
        notes_in_dict = sp.parseNotes(rxn)
        notes_in_dict["GENE_ASSOCIATION"] = [study_ga]
        notes = "<body xmlns=\"http://www.w3.org/1999/xhtml\">"
        for k,v_list in notes_in_dict.items():
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
    config.set('DATABASE_PATHS', 'database_ref_path', '/home/data/database/BIOCYC/Metacyc/22.0_enhanced/metacyc_22.0_enhanced.padmet')
    config.set('DATABASE_PATHS', 'mnx_rxn_path', '/data/database/MNX/reac_xref.csv')
    config.set('DATABASE_PATHS', 'mnx_cpd_path', '/data/database/MNX/chem_xref.csv')
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
    config.set('TOOL_PATHS', 'orthofinder_bin_path', '/programs/OrthoFinder-2.2.6/orthofinder')
    config.set('TOOL_PATHS', 'padmet_utils_path', '/programs/padmet-utils')
    config.add_section('VAR')
    config.set('VAR', 'study_from_annot_prefix', 'output_pathwaytools_')

    # Writing our configuration file to 'example.cfg'
    with open(config_file_path, 'wb') as configfile:
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

if __name__ == "__main__":
    main()