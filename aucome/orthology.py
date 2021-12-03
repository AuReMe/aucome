#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
usage:
     aucome orthology --run=ID [--sequence_search_prg=STR] [--cpu=INT] [-v] [--vv] [--filtering] [--threshold=FLOAT] [--union] [--intersection]
     
options:
     --run=ID    Pathname to the comparison workspace.
     --sequence_search_prg=STR    Sequence search program for Orthofinder [Default: diamond].
         Options: blast, mmseqs, blast_gz, diamond, and diamond_ultra_sens.
     --cpu=INT     Number of cpu to use for the multiprocessing (if none use 1 cpu).
     -v     Verbose.
     --vv    Very verbose.
     --filtering     Use a filter to limit propagation, by default it is 0.05, if you want to modify the value use --threshold.
     --threshold=FLOAT     Threshold of the filter to limit propagation to use with the --filtering argument.
    --union          Use the union filter between five threshold values [0.01, 0.05, 0.1, 0.15, 0.2] to limit propagation, to use with the --filtering argument.
    --intersection   Use the intersection filter between five threshold values [0.01, 0.05, 0.1, 0.15, 0.2] to limit propagation, to use with the --filtering argument.
"""

import csv
import docopt
import os
import subprocess
import shutil
import sys
import time

from padmet.classes import PadmetSpec
from padmet.utils.exploration import convert_sbml_db
from padmet.utils.connection import extract_orthofinder
from padmet.utils.connection import sbml_to_padmet, sbmlGenerator

from aucome.utils import parse_config_file
from multiprocessing import Pool


def command_help():
    print(docopt.docopt(__doc__))


def orthology_parse_args(command_args):
    args = docopt.docopt(__doc__, argv=command_args)
    run_id = args['--run']
    sequence_search_prg = args['--sequence_search_prg']
    cpu = args['--cpu']
    verbose = args['-v']
    veryverbose = args['--vv']
    filtering = args['--filtering']
    threshold = args['--threshold']
    union = args['--union']
    intersection = args['--intersection']
    filtering_threshold_list = []
    
    if filtering:
        if threshold:
            try:
                filtering_threshold_list.append(float(threshold))
            except:
                sys.exit('filtering_threshold value must be a float between 0 and 1.')       
        else:
            filtering_threshold_list.append(0.05)      
        if union or intersection:
            filtering_threshold_list = [0.01, 0.05, 0.1, 0.15, 0.2]
    else:
        if threshold:
            sys.exit('--threshold must be used with --filtering.')
        if union:
            sys.exit('--union must be used with --filtering.')
        if intersection:
            sys.exit('--intersection must be used with --filtering.')
        
    if cpu:
        nb_cpu_to_use = int(cpu)
    else:
        nb_cpu_to_use = 1

    if veryverbose and not verbose:
        verbose = veryverbose

    run_orthology(run_id, sequence_search_prg, nb_cpu_to_use, filtering_threshold_list, union, intersection, verbose, veryverbose)


def run_orthology(run_id, sequence_search_prg, nb_cpu_to_use, filtering_threshold_list, union, intersection, verbose, veryverbose=None):
    print('--- Running orthology step ---')
    orthology_start_time = time.time()
    aucome_pool = Pool(nb_cpu_to_use)
    config_data = parse_config_file(run_id)

    orthofinder_wd_path = config_data['orthofinder_wd_path']
    orthofinder_bin_path = config_data['orthofinder_bin_path']
    orthofinder_sbml_path = config_data['orthofinder_sbml_path']
    orthofinder_padmet_path = config_data['orthofinder_padmet_path']
    orthofinder_filtered_path = config_data['orthofinder_filtered_path']
    studied_organisms_path = config_data['studied_organisms_path']
    padmet_from_annotation_path = config_data['padmet_from_annotation_path']
    database_path = config_data['database_path']

    all_study_name = set(next(os.walk(studied_organisms_path))[1])
    
    all_study_faa = dict([(study_name, "{0}/{1}/{1}.faa".format(studied_organisms_path, study_name))
                          if os.path.isfile("{0}/{1}/{1}.faa".format(studied_organisms_path, study_name))
                          else (study_name, '')
                          for study_name in all_study_name])

    #check if Orthofinder already run, if yes, get the last workdir.
    try:
    #    if orthogroups:
    #        orthodata_path = max(["%s/%s" %(x[0], 'Orthogroups/Orthogroups.tsv') for x in os.walk(orthofinder_wd_path) if 'Orthogroups' in x[1]])
    #    else:
        orthodata_path = max(["%s/%s" %(x[0], 'Orthologues') for x in os.walk(orthofinder_wd_path) if 'Orthologues' in x[1]])
    except ValueError:
        if verbose:
            print("Unable to find file Orthogroups.tsv in {0}, need to run Orthofinder...".format(orthofinder_wd_path))
        orthodata_path = None

    # If there is already an orthofinder result, check if all the species are in it.
    if orthodata_path:
        wd_orthodata_path = max(["%s/%s" %(x[0], 'WorkingDirectory') for x in os.walk(orthofinder_wd_path) if 'WorkingDirectory' in x[1]])

        input_fasta = [fasta_name for fasta_name in all_study_faa]
        already_analysed_fasta = [fasta_name.replace('.faa', '') for fasta_name in os.listdir(orthofinder_wd_path) if fasta_name != 'OrthoFinder']

        # If there is missing species, rerun OrthoFinder to add the missing species.
        if len(already_analysed_fasta) != len(input_fasta):
            fasta_to_adds = set(input_fasta) - set(already_analysed_fasta)
            tmp_folder = orthofinder_wd_path + '/tmp/'
            os.mkdir(tmp_folder)

            print('There is missing species in orthofinder: ' + ','.join(list(fasta_to_adds)))
            print('Rerun OrthoFinder on them using the old results from '+ orthofinder_wd_path)

            for name, faa_path in list(all_study_faa.items()):
                if name in fasta_to_adds:
                    if not os.path.isfile("{0}/{1}.faa".format(tmp_folder, name)):
                        if verbose:
                            print("Copying {0}'s faa to {1}".format(name, tmp_folder))
                    if faa_path != "" and os.path.isfile(faa_path):
                        shutil.copyfile(faa_path, tmp_folder + '/' + name + '.faa')
                    else:
                        sys.exit("Missing fasta for " + name + ", use 'aucome check' to input fasta from genbank.")

            orthofinder_result_path = orthofinder_wd_path + '/OrthoFinder/'
            start_time = time.time()
            cmds = [orthofinder_bin_path, "-b", wd_orthodata_path, "-f", tmp_folder,
                    "-t", str(nb_cpu_to_use), "-S", sequence_search_prg]
                       
            subprocess.call(cmds)
            end_time = (time.time() - start_time)
            integer_part, decimal_part = str(end_time).split('.')
            end_time = ".".join([integer_part, decimal_part[:3]])
            if verbose:
                print("Orthofinder done in: %ss" %end_time)

            for fasta_name in os.listdir(tmp_folder):
                shutil.copyfile(tmp_folder + '/' + fasta_name, orthofinder_wd_path + '/' + fasta_name)

            shutil.rmtree(tmp_folder)

            # Replace the old OrthoFinder results folder with the new one.
            new_orthofidner_path = max(["%s/%s" %(x[0], 'OrthoFinder') for x in os.walk(wd_orthodata_path) if 'OrthoFinder' in x[1]])

            orthofinder_tmp = orthofinder_wd_path + '/OrthoFinder_tmp'

            shutil.copytree(new_orthofidner_path, orthofinder_tmp)
            shutil.rmtree(orthofinder_result_path)
            shutil.copytree(orthofinder_tmp, orthofinder_result_path)
            shutil.rmtree(orthofinder_tmp)

            #if orthogroups:
            #    orthodata_path = max(["%s/%s" %(x[0], 'Orthogroups/Orthogroups.tsv') for x in os.walk(orthofinder_wd_path) if 'Orthogroups' in x[1]])
            #else:
            orthodata_path = max(["%s/%s" %(x[0], 'Orthologues') for x in os.walk(orthofinder_wd_path) if 'Orthologues' in x[1]])

            # Clean sbml/padmet/padmet filtered to recreate them with the new data.
            for sbml_folder in os.listdir(orthofinder_sbml_path):
                shutil.rmtree(orthofinder_sbml_path + '/' + sbml_folder)

            for padmet_folder in os.listdir(orthofinder_padmet_path):
                os.remove(orthofinder_padmet_path + '/' + padmet_folder)

            for filtered_padmet_folder in os.listdir(orthofinder_filtered_path):
                os.remove(orthofinder_filtered_path + '/' + filtered_padmet_folder)

    # If there is no OrthoFinder results folder run OrthoFinder on all the fasta.
    elif not orthodata_path:
        for name, faa_path in list(all_study_faa.items()):
            if not os.path.isfile("{0}/{1}.faa".format(orthofinder_wd_path, name)):
                if verbose:
                    print("Copying {0}'s faa to {1}".format(name, orthofinder_wd_path))
            if faa_path != "" and os.path.isfile(faa_path):
                shutil.copyfile(faa_path, orthofinder_wd_path + '/' + name + '.faa')
            else:
                sys.exit("Missing fasta for " + name + ", use 'aucome check' to input fasta from genbank.")
                
        if verbose:
            print("Running Orthofinder on %s cpu" %nb_cpu_to_use)

        start_time = time.time()
        cmds = [orthofinder_bin_path, "-f", orthofinder_wd_path, "-t", str(nb_cpu_to_use), "-S", sequence_search_prg]
        subprocess.call(cmds)
        end_time = (time.time() - start_time)
        integer_part, decimal_part = str(end_time).split('.')
        end_time = ".".join([integer_part, decimal_part[:3]])
        if verbose:
            print("Orthofinder done in: %ss" %end_time)

        orthofinder_result_path = orthofinder_wd_path + '/OrthoFinder/'
        if os.path.exists(orthofinder_result_path):
            orthologues_result_path = orthofinder_result_path + os.listdir(orthofinder_result_path)[0] + '/Orthologues'
            if os.path.exists(orthologues_result_path):
                orthodata_path = max(["%s/%s" %(x[0], 'Orthologues') for x in os.walk(orthofinder_wd_path) if 'Orthologues' in x[1]])
            else:
                sys.exit('There was an error with OrthoFinder, there is no results in ' + orthologues_result_path)
        else:
            sys.exit('Missing OrthoFinder folder in ' + orthofinder_wd_path)

    if verbose:
        print("Parsing Orthofinder output %s" %orthodata_path)

    if verbose:
        print("Start sbml creation...")
    all_dict_data = []
    for study_name in all_study_name:
        output_sbml = os.path.join(orthofinder_sbml_path, study_name)
        if os.path.exists(output_sbml):
            print(output_sbml + " already exists, delete it if you want to relaunch ortholog creation.")
        else:
            dict_data = {'sbml': run_id, 'orthodata_path': orthodata_path,
                         'study_name': study_name, 'verbose': verbose,
                         'veryverbose': veryverbose, 'output': output_sbml}
            all_dict_data.append(dict_data)

    start_time = time.time()
    aucome_pool.map(orthogroup_to_sbml, all_dict_data)
    end_time = (time.time() - start_time)
    integer_part, decimal_part = str(end_time).split('.')
    end_time = ".".join([integer_part, decimal_part[:3]])
    if verbose:
        print("Orthofinder output parsed in: %ss" %end_time)
    """
    #check database, mapping to metacyc ???
    data_convert_sbml_db = []
    for dict_data in all_dict_data:
        tmp_dict_data = {'sbml': orthology_based_path + '/' + study_name, 
                         'mnx_rxn_path': mnx_rxn_path, 'mnx_cpd_path': mnx_cpd_path, 'verbose': verbose}
        data_convert_sbml_db.append(tmp_dict_data)
        
    aucome_pool.map(_convert_sbml_db, data_convert_sbml_db)
    """

    if verbose:
        if len(filtering_threshold_list)>0 :
            print("Start padmet creation and filtering...")
        else:
            print("Start padmet creation...")

    multiprocessing_datas = []
    update_padmet_datas = []
    for sbml in os.listdir(orthofinder_sbml_path):
        input_pwt_padmet = padmet_from_annotation_path + '/output_pathwaytools_' + sbml + '.padmet'
        output_padmet = orthofinder_padmet_path + '/' + sbml + '.padmet'
        if os.path.exists(output_padmet):
            print(output_padmet + " already exists, delete it if you want to relaunch ortholog creation.")
        else:
            multiprocessing_datas.append([sbml, orthofinder_sbml_path, input_pwt_padmet,
                            database_path, output_padmet, orthodata_path,
                            orthofinder_filtered_path, filtering_threshold_list, verbose, veryverbose])
            update_padmet_datas.append([sbml, orthodata_path, output_padmet, verbose])

    start_time = time.time()
    aucome_pool.starmap(orthology_to_padmet, multiprocessing_datas)

    # Add the orthologs to the padmets.
    if verbose:
        print("Updating padmets...")
    aucome_pool.starmap(addOrthologyInPadmet, update_padmet_datas)

    if len(filtering_threshold_list)>0:
        filter_propagation(orthofinder_padmet_path, orthofinder_filtered_path, aucome_pool, filtering_threshold_list, union, intersection, verbose)
   
    end_time = (time.time() - start_time)
    integer_part, decimal_part = str(end_time).split('.')
    end_time = ".".join([integer_part, decimal_part[:3]])
    if verbose:
        print("Padmet created in: %ss" %end_time)

    aucome_pool.close()
    aucome_pool.join()

    orthology_end_time = (time.time() - orthology_start_time)
    integer_part, decimal_part = str(orthology_end_time).split('.')
    orthology_time = ".".join([integer_part, decimal_part[:3]])

    if verbose:
        print("--- orthology step done in: %ss ---" %orthology_time)


def _convert_sbml_db(data_convert_sbml_db):
    
    sbml_file = data_convert_sbml_db['sbml']
    verbose = data_convert_sbml_db['verbose']
    mnx_rxn_path = data_convert_sbml_db['mnx_rxn_path']
    mnx_cpd_path = data_convert_sbml_db['mnx_cpd_path']

    if os.path.isfile(sbml_file):
        dict_file = "{0}_dict.tsv".format(os.path.splitext(sbml_file)[0])
        if not os.path.exists(dict_file):
            db_ref = convert_sbml_db.check_sbml_db(sbml_file, "reaction", mnx_reac_file=mnx_rxn_path, verbose=True)[0]
            if verbose:
                print("%s: %s" %(os.path.basename(sbml_file), db_ref))
            if db_ref.lower() != "metacyc":
                if verbose:
                    print("Creating id mapping file: %s" %dict_file)
                convert_sbml_db.map_sbml(sbml_file, "reaction", "metacyc", dict_file, verbose=verbose, mnx_reac_file=mnx_rxn_path, mnx_chem_file=mnx_cpd_path)


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
    veryverbose = dict_data['veryverbose']
    #orthogroups = dict_data['orthogroups']

    if verbose:
        print('Create sbml for ' + study_name)

    if veryverbose:
        extract_orthofinder_verbose = veryverbose
    else:
        extract_orthofinder_verbose = False

    all_model_sbml = extract_orthofinder.get_sbml_files(sbml, workflow="aucome", verbose=extract_orthofinder_verbose)
    #if orthogroups:
    #    extract_orthofinder.orthogroups_to_sbml(orthodata_path, all_model_sbml, output, study_name, extract_orthofinder_verbose)
    #else:
    extract_orthofinder.orthologue_to_sbml(orthodata_path, all_model_sbml, output, study_name, extract_orthofinder_verbose)


def orthology_to_padmet(sbml, orthofinder_sbml_path, input_pwt_padmet,
                         database_path, output_padmet, orthodata_path,
                         orthofinder_filtered_path, filtering_threshold_list,
                        verbose, veryverbose):
    source_tool = "ORTHOFINDER"
    source_category = "ORTHOLOGY"

    if verbose:
        print('Create padmet from sbml for ' + sbml)
    if veryverbose:
        sbml_padmet_verbose = veryverbose
    else:
        sbml_padmet_verbose = False

    sbml_to_padmet.sbml_to_padmetSpec(orthofinder_sbml_path + '/' + sbml,
                                      input_pwt_padmet,
                                      padmetRef_file=database_path,
                                      output=output_padmet,
                                      source_tool=source_tool,
                                      source_category=source_category,
                                      verbose=sbml_padmet_verbose)


def addOrthologyInPadmet(study_id, orthodata_path, output_padmet, verbose):
    """
    Add orthologs information to a padmet file.
    Args:
        study_id (str): ID of the species which padmet will be analyzed
        orthodata_path (str): path to Orthologues files
        output_padmet (str): path to the output padmet file
        verbose (boolean): verbose
    """
    # Read orthologs files and create a dictionary containing orthologs
    # infomations.
    all_orgs = set()
    all_orthologue_files = []
    # Process of orthologues.
    for _path, _folders, _files in os.walk(orthodata_path):
        for _folder in _folders:
            all_orgs.add(_folder.replace("Orthologues_","").lower())
        for _file in _files:
            _filename = os.path.splitext(_file)[0]
            if _filename.startswith(study_id):
                all_orthologue_files.append(os.path.join(_path,_file))
        dict_orthologues = {}

    for org in all_orgs:
        dict_orthologues[org] = dict()

    #dict_orthologue: k = org_id, v = dict: k = gene_id, v = dict: k = org_id, v = set of gene orthologue id
    for orthologue_file in [i for i in all_orthologue_files if "__v__" in i]:
        with open(orthologue_file, 'r') as csvfile:
            reader = csv.DictReader(csvfile, delimiter = "\t")
            orgs = list(reader.fieldnames)
            orgs.remove('Orthogroup')
            org_A, org_B = orgs
            for row in reader:
                gene_ids_A = [gene_id.split("_isoform")[0] for gene_id in row[org_A].split(", ")]
                gene_ids_B = [gene_id.split("_isoform")[0] for gene_id in row[org_B].split(", ")]
                for gene_id_A in gene_ids_A:
                    if gene_id_A not in dict_orthologues[org_A.lower()].keys():
                        dict_orthologues[org_A.lower()][gene_id_A] = dict()
                    dict_orthologues[org_A.lower()][gene_id_A][org_B.lower()] = set(gene_ids_B)

    padmet_file = os.path.basename(output_padmet)
    if verbose:
        print("%s..."%padmet_file)
    org_id = os.path.splitext(padmet_file)[0]

    org_id = org_id.lower()
    padmet = PadmetSpec(output_padmet)
    for linked_rlt in [rlt for rlt in padmet.getAllRelation() if rlt.type == "is_linked_to"]:
        gene_id = linked_rlt.id_out
        for index, src in enumerate(linked_rlt.misc["SOURCE:ASSIGNMENT"]):
            if 'GENOME:' not in src:
                ortho_org_id = src.replace("OUTPUT_ORTHOFINDER_FROM_","").split(':')[0]
                ortho_org_id_low = ortho_org_id.lower()
                org_id_orthologues = dict_orthologues[org_id][gene_id]
                ortho_genes_id = ";".join(org_id_orthologues[ortho_org_id_low])
                new_src = "OUTPUT_ORTHOFINDER_FROM_%s:%s"%(ortho_org_id, ortho_genes_id)
                linked_rlt.misc["SOURCE:ASSIGNMENT"][index] = new_src
    padmet.generateFile(output_padmet)


def filter_propagation(padmet_folder, output_folder, aucome_pool, filtering_threshold_list, union=None, intersection=None, verbose=None):
    propagation_to_remove_file = os.path.join(output_folder, "propagation_to_remove.tsv")
    reactions_to_remove_file = os.path.join(output_folder, 'reactions_to_remove.tsv')
    padmet_output_folder = output_folder

    if verbose:
        print("Extracting all the relations gene-reaction...")
    dict_rxn_orgs_genes, dict_rxn_ec = extractRGL(padmet_folder, aucome_pool)
    if verbose:
        print("Extracting all the gene propagations...")
    dict_rxn_org_gene_propagation = extractPropagation(dict_rxn_orgs_genes)
    if verbose:
        print("Writing the file propagation_to_remove...")
    dict_rxn_org_gene_propag_to_remove = extractPropagationToRemove(dict_rxn_org_gene_propagation, output=propagation_to_remove_file, orthology_threshold_list=filtering_threshold_list, union=union, intersection=intersection)
    if verbose:
        print("Cleaning the Padmet files and writing the reactions_to_remove_file file...")
    cleanPadmet(dict_rxn_org_gene_propag_to_remove, dict_rxn_ec, padmet_folder,
                padmet_output_folder, reactions_to_remove_file, aucome_pool)


def merge_result_dict(first_dict, second_dict):
    for rxn_id in second_dict:
        if rxn_id not in first_dict:
            first_dict[rxn_id] = dict()
        for org_id in second_dict[rxn_id]:
            if org_id not in first_dict[rxn_id]:
                first_dict[rxn_id][org_id] = dict()
            for gene_id in second_dict[rxn_id][org_id]:
                if gene_id not in first_dict[rxn_id][org_id]:
                    first_dict[rxn_id][org_id][gene_id] = second_dict[rxn_id][org_id][gene_id]

def extractRGL(padmet_folder, aucome_pool):
    """
    extract reactions genes relations.
    It reads all Padmet files in padmet_folder, then it creates three
    dictionaries: dict_rnx_orgs_genes, dict_rnx_ec, dict_org_genes.
    dict_rxn_org_gene: reaction & organism & assigned genes & orthologuous genes.
    dict_rxn_ec[rxn_id] = ec, simple dictionary
    return dict {rxn_id:{org_id:{'FROM-PTOOL': bool, gene_id:set of sources (ex ortho_org_id:ortho_genes)}}}
    """
    dict_rxn_orgs_genes = {}
    dict_rxn_ec = {}
    multiprocessing_datas = []
    for padmet_file in next(os.walk(padmet_folder))[2]:
        multiprocessing_datas.append([padmet_file, padmet_folder])

    multiprocessing_results = aucome_pool.starmap(mp_extractRGL, multiprocessing_datas)

    for multiprocessing_result in multiprocessing_results:
        merge_result_dict(dict_rxn_orgs_genes, multiprocessing_result[0])
        dict_rxn_ec.update(multiprocessing_result[1])

    return dict_rxn_orgs_genes, dict_rxn_ec


def mp_extractRGL(padmet_file, padmet_folder):
    padmet_path = os.path.join(padmet_folder, padmet_file)
    org_id = os.path.splitext(padmet_file)[0].upper()
    padmet = PadmetSpec(padmet_path)
    dict_rxn_orgs_genes = {}
    dict_rxn_ec = {}
    for node in padmet.dicOfNode.values():
        if node.type == "reaction":
            rxn_id = node.id
            if rxn_id not in dict_rxn_orgs_genes:
                dict_rxn_orgs_genes[rxn_id] = dict()

            if org_id not in dict_rxn_orgs_genes[rxn_id]:
                dict_rxn_orgs_genes[rxn_id][org_id] = dict()
            if 'EC-NUMBER' in node.misc:
                dict_rxn_ec[rxn_id] = node.misc['EC-NUMBER'][0]
            #rxn_id = "PHOSGLYPHOS-RXN"
            from_ptool = False
            for is_linked_rlt in [rlt for rlt in padmet.dicOfRelationIn[rxn_id] if rlt.type == "is_linked_to"]:
                gene_id = is_linked_rlt.id_out
                all_sources = {src.replace("OUTPUT_ORTHOFINDER_FROM_","")  for src in is_linked_rlt.misc["SOURCE:ASSIGNMENT"]}
                if any (src for src in all_sources if src.startswith("GENOME")):
                    from_ptool = True
                dict_rxn_orgs_genes[rxn_id][org_id][gene_id] = all_sources
            dict_rxn_orgs_genes[rxn_id][org_id]["FROM-PTOOL"] = from_ptool

    return dict_rxn_orgs_genes, dict_rxn_ec


def extractPropagation(dict_rxn_orgs_genes):
    """
    Propagations are extracted from Padmet networks. Then propagations
    are split bewteen two groups. Those who are orthologues to gene-reaction
    associations that come from Pathway Tools, and those whose associations
    do not come from Pathway Tools.
    """
    dict_rxn_org_gene_propagation = {}
    for rxn_id, rxn_data in dict_rxn_orgs_genes.items():
        dict_rxn_org_gene_propagation[rxn_id] = dict()
        for org_id, org_data in rxn_data.items():
            if org_id not in dict_rxn_org_gene_propagation[rxn_id].keys():
                dict_rxn_org_gene_propagation[rxn_id][org_id] = dict()
            #because of genes without gene link sources
            if len(org_data.keys()) > 1:
                for gene_id, gene_data in org_data.items():
                    if gene_id != "FROM-PTOOL":
                        gene_id = gene_id
                        if any (src for src in gene_data if src.startswith("GENOME")):
                            is_from_ptool = True
                        else:
                            is_from_ptool = False
                        for src in gene_data:
                            if not src.startswith("GENOME:"):
                                ortho_org_id, ortho_genes_ids = src.split(":")
                                ortho_genes_ids = ortho_genes_ids.split(";")
                                if ortho_org_id not in dict_rxn_org_gene_propagation[rxn_id].keys():
                                    dict_rxn_org_gene_propagation[rxn_id][ortho_org_id] = dict()
                                for ortho_gene_id in ortho_genes_ids:
                                    if gene_id in dict_rxn_org_gene_propagation[rxn_id][org_id]:
                                        dict_rxn_org_gene_propagation[rxn_id][org_id][gene_id]["propagation_from_ptool"].add((ortho_org_id, ortho_gene_id))
                                    else:
                                        dict_rxn_org_gene_propagation[rxn_id][org_id][gene_id] = {"propagation_to_ptool":set(),"propagation_to_not_ptool":set(), "propagation_from_ptool": set([(ortho_org_id, ortho_gene_id)])}

                                    if ortho_gene_id not in dict_rxn_org_gene_propagation[rxn_id][ortho_org_id].keys():
                                        dict_rxn_org_gene_propagation[rxn_id][ortho_org_id][ortho_gene_id] = {"propagation_to_ptool":set(),"propagation_to_not_ptool":set(), "propagation_from_ptool": set()}
                                    if is_from_ptool:
                                        dict_rxn_org_gene_propagation[rxn_id][ortho_org_id][ortho_gene_id]["propagation_to_ptool"].add((org_id, gene_id))
                                    else:
                                        dict_rxn_org_gene_propagation[rxn_id][ortho_org_id][ortho_gene_id]["propagation_to_not_ptool"].add((org_id, gene_id))

    return dict_rxn_org_gene_propagation


def extractPropagationToRemove(dict_rxn_org_gene_propagation, output,
                               ptool_threshold=0,
                               orthology_threshold_list=[0.05], union=None,
                               intersection=None):
    """
    Using ptool_threshold and orthology_threshold, this function select the 
    propagations to remove. These propagation are written in 
    propagation_to_remove.tsv.
    """
    header = ["reaction_id", "org_id", "gene_id"]
    dict_rxn_org_gene_propag_to_remove = dict()
    maximum = 5
    with open(output, 'w') as csvfile:
        writer = csv.DictWriter(csvfile, header, delimiter="\t")
        writer.writeheader()
        for rxn_id ,rxn_data in dict_rxn_org_gene_propagation.items():
            nb_org_prop = len(rxn_data.keys())-1
            if nb_org_prop:
                for org_id, org_data in rxn_data.items():
                    gene_ids_to_remove = set()
                    gene_id_tmp = ""
                    for gene_id, gene_data in org_data.items():
                        count = 0
                        # At this moment filter is as 20/N with 0.05
                        for orthology_threshold in orthology_threshold_list:
                            inverse_orthology_threshold = 1/orthology_threshold
                            not_ptool_threshold = round(max(inverse_orthology_threshold/nb_org_prop, orthology_threshold*nb_org_prop),0)
                            if len({i[0] for i in gene_data["propagation_to_ptool"]}) <= ptool_threshold:
                                if len({i[0] for i in gene_data["propagation_to_not_ptool"]}) >= not_ptool_threshold:                                    
                                    if intersection:
                                        count+=1 
                                        gene_id_tmp = gene_id
                                    elif not gene_id in gene_ids_to_remove:
                                        gene_ids_to_remove.add(gene_id)
                                        dict_rxn_org_gene_propag_to_remove = remove_gene(gene_ids_to_remove, dict_rxn_org_gene_propag_to_remove, rxn_id, org_id, writer)
                        if intersection and count == maximum:
                            gene_ids_to_remove.add(gene_id)
                            dict_rxn_org_gene_propag_to_remove = remove_gene(gene_ids_to_remove, dict_rxn_org_gene_propag_to_remove, rxn_id, org_id, writer)
    return dict_rxn_org_gene_propag_to_remove


# This function is called in extractPropagationToRemove().
def remove_gene(gene_ids_to_remove, dict_rxn_org_gene_propag_to_remove, rxn_id,
                org_id, writer):
    if gene_ids_to_remove:
        genes_ids = ";".join(gene_ids_to_remove)
        if rxn_id in dict_rxn_org_gene_propag_to_remove:
            dict_rxn_org_gene_propag_to_remove[rxn_id][org_id] = gene_ids_to_remove
        else:
            dict_rxn_org_gene_propag_to_remove[rxn_id] = {org_id: gene_ids_to_remove}
            line = {"reaction_id": rxn_id, "org_id": org_id, "gene_id": genes_ids}
            writer.writerow(line)
    return dict_rxn_org_gene_propag_to_remove


def merge_dict_org_rxn_clean(dict_org_rxn_clean, tmp_dict_org_rxn_clean,
                             dict_rxn_org_gene_propag_to_remove):
    for org_id in tmp_dict_org_rxn_clean:
        if org_id not in dict_org_rxn_clean:
            dict_org_rxn_clean[org_id] = dict()
        for rxn_id in tmp_dict_org_rxn_clean[org_id]:
            if rxn_id in dict_rxn_org_gene_propag_to_remove.keys():
                if rxn_id not in dict_org_rxn_clean[org_id]:
                    dict_org_rxn_clean[org_id][rxn_id] = dict()
                for gene_id in tmp_dict_org_rxn_clean[org_id][rxn_id]:
                    if gene_id not in dict_org_rxn_clean[org_id][rxn_id]:
                        dict_org_rxn_clean[org_id][rxn_id][gene_id] = tmp_dict_org_rxn_clean[org_id][rxn_id][gene_id]


def create_dict_org_rxn_clean(padmet_file, padmet_folder,
                              dict_rxn_org_gene_propag_to_remove):
    padmet_path = os.path.join(padmet_folder, padmet_file)
    org_id = os.path.splitext(padmet_file)[0].upper()
    padmet = PadmetSpec(padmet_path)
    dict_org_rxn_clean = dict()
    dict_org_rxn_clean[org_id] = dict()
    for rxn_id in [node.id for node in padmet.dicOfNode.values() if node.type == "reaction" and node.id in dict_rxn_org_gene_propag_to_remove.keys()]:
        dict_org_rxn_clean[org_id][rxn_id] = dict()
        for is_linked_rlt in [rlt for rlt in padmet.dicOfRelationIn[rxn_id] if rlt.type == "is_linked_to"]:
            gene_id = is_linked_rlt.id_out
            all_sources = {src.replace("OUTPUT_ORTHOFINDER_FROM_","")  for src in is_linked_rlt.misc["SOURCE:ASSIGNMENT"]}
            new_sources = set()
            for src in all_sources:
                if src.startswith("GENOME:"):
                    new_sources.add(src)
                else:
                    ortho_org_id, ortho_genes_ids = src.split(":")
                    ortho_genes_ids = set(ortho_genes_ids.split(";"))
                    if ortho_org_id not in dict_rxn_org_gene_propag_to_remove[rxn_id].keys():
                        new_sources.add(src)
                    else:
                        new_ortho_genes_ids = set()
                        for ortho_gene_id in ortho_genes_ids:
                            if ortho_gene_id not in dict_rxn_org_gene_propag_to_remove[rxn_id][ortho_org_id]:
                                new_ortho_genes_ids.add(ortho_gene_id)
                        if new_ortho_genes_ids:
                            new_src = "%s:%s"%(ortho_org_id, ";".join(new_ortho_genes_ids))
                            new_sources.add(new_src)
            dict_org_rxn_clean[org_id][rxn_id][gene_id] = new_sources

    return dict_org_rxn_clean


# This function is called in the cleanPadmet() function. 
def delete_propagation(padmet_file, padmet_folder, dict_org_rxn_clean,
                       output_folder):
    padmet_path = os.path.join(padmet_folder, padmet_file)
    org_id = os.path.splitext(padmet_file)[0].upper()
    padmet = PadmetSpec(padmet_path)
    output = os.path.join(output_folder, padmet_file)
    nb_rxn_removed = 0
    for rxn_id, rxn_data in dict_org_rxn_clean[org_id].items():
        if any(rxn_data.keys()):
            if not any(rxn_data.values()):
                nb_rxn_removed += 1
                padmet.delNode(rxn_id)
            else:
                for gene_id, gene_data in rxn_data.items():
                    is_linked_rlt = [rlt for rlt in padmet.dicOfRelationOut[gene_id] if rlt.id_in == rxn_id and rlt.id_out == gene_id][0]
                    if not gene_data:
                        #remove relation
                        padmet._delRelation(is_linked_rlt)
                    else:
                        #update relation.misc[src:assgn], /!\ MAJ et source edition
                        is_linked_rlt.misc["SOURCE:ASSIGNMENT"] = ["OUTPUT_ORTHOFINDER_FROM_%s"%src for src in gene_data]
    print("Removing %s in %s"%(nb_rxn_removed, org_id))
    padmet.generateFile(output)


def cleanPadmet(dict_rxn_org_gene_propag_to_remove, dict_rxn_ec, padmet_folder,
                output_folder, reactions_to_remove_file, aucome_pool):
    """
    It cleans the Padmet files and it writes the reactions_to_remove_file file.
    """
    dict_org_rxn_clean = dict()
    create_multiprocessing_datas = []
    for padmet_file in next(os.walk(padmet_folder))[2]:
        create_multiprocessing_datas.append([padmet_file, padmet_folder, dict_rxn_org_gene_propag_to_remove])

    multiprocessing_results = aucome_pool.starmap(create_dict_org_rxn_clean, create_multiprocessing_datas)

    for multiprocessing_result in multiprocessing_results:
        merge_dict_org_rxn_clean(dict_org_rxn_clean, multiprocessing_result, dict_rxn_org_gene_propag_to_remove)

    remove_multiprocessing_datas = []
    for padmet_file in next(os.walk(padmet_folder))[2]:
        remove_multiprocessing_datas.append([padmet_file, padmet_folder, dict_org_rxn_clean, output_folder])

    aucome_pool.starmap(delete_propagation, remove_multiprocessing_datas)

    with open(reactions_to_remove_file, 'w') as csvfile:
        header = ["org_id","reaction_id", "ec-number", "gene_id"]
        writer = csv.DictWriter(csvfile, header, delimiter="\t")
        writer.writeheader()
        list_count = []
        for org_id, org_data in dict_org_rxn_clean.items():
            to_remove = list()
            for rxn_id, rxn_data in org_data.items():
                if any(rxn_data.keys()) and not any(rxn_data.values()):
                    to_remove.append(rxn_id)
                    if rxn_id in dict_rxn_ec:
                        ec = dict_rxn_ec[rxn_id]
                    else:
                        ec = ''
                    gene_id = dict_org_rxn_clean[org_id][rxn_id]
                    line = {"org_id":org_id, "reaction_id":rxn_id,
                            "ec-number":ec, "gene_id":";".join(gene_id)}
                    writer.writerow(line)
                    list_count.append(rxn_id)
    return dict_org_rxn_clean
