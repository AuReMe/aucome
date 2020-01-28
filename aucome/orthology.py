#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
usage:
    aucome orthology --run=ID [-S=STR] [--orthogroups] [--cpu=INT] [-v] [--filtering]

options:
    --run=ID    Pathname to the comparison workspace.
    --orthogroups    Use Orthogroups instead of Orthologues after Orthofinder.
    -S=STR    Sequence search program for Orthofinder [Default: diamond].
        Options: blast, mmseqs, blast_gz, diamond
    --cpu=INT     Number of cpu to use for the multiprocessing (if none use 1 cpu).
    -v     Verbose.
    --filtering     Use a filter to limit propagation.
"""

import csv
import docopt
import os
import subprocess
import shutil
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
    orthogroups = args['--orthogroups']
    sequence_search_prg = args['-S']
    verbose = args['-v']
    filtering = args['--filtering']

    if args["--cpu"]:
        nb_cpu_to_use = int(args["--cpu"])
    else:
        nb_cpu_to_use = 1

    run_orthology(run_id, orthogroups, sequence_search_prg, nb_cpu_to_use, filtering, verbose)

def run_orthology(run_id, orthogroups, sequence_search_prg, nb_cpu_to_use, filtering, verbose):
    aucome_pool = Pool(nb_cpu_to_use)

    config_data = parse_config_file(run_id)

    orthofinder_wd_path = config_data['orthofinder_wd_path']
    orthofinder_bin_path = config_data['orthofinder_bin_path']
    orthofinder_sbml_path = config_data['orthofinder_sbml_path']
    orthofinder_padmet_path = config_data['orthofinder_padmet_path']
    orthofinder_filtered_path = config_data['orthofinder_filtered_path']
    studied_organisms_path = config_data['studied_organisms_path']
    model_organisms_path = config_data['model_organisms_path']
    padmet_from_annotation_path = config_data['padmet_from_annotation_path']
    database_path = config_data['database_path']

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
                shutil.copyfile(faa_path, orthofinder_wd_path + '/' + name + '.faa')

        for name, faa_path in list(all_model_faa.items()):
            if not os.path.isfile("{0}/{1}.faa".format(orthofinder_wd_path, name)):
                if verbose:
                    print("Copying {0}'s faa to {1}".format(name, orthofinder_wd_path))
                shutil.copyfile(faa_path, orthofinder_wd_path + '/' + name + '.faa')

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
                    'output': orthofinder_sbml_path + '/' + study_name}
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
    """
    aucome_pool.close()
    aucome_pool.join()

    source_tool = "ORTHOFINDER"
    source_category = "ORTHOLOGY"
    for sbml in os.listdir(orthofinder_sbml_path):
        sbml_to_padmet.sbml_to_padmetSpec(orthofinder_sbml_path + '/' + sbml,
                                        padmet_from_annotation_path + '/output_pathwaytools_' + sbml + '.padmet',
                                        padmetRef_file=database_path,
                                        output=orthofinder_padmet_path + '/' + sbml + '.padmet',
                                        source_tool=source_tool, source_category=source_category, verbose=verbose)

    addOrthologyInPadmet(orthodata_path, orthofinder_padmet_path, orthofinder_padmet_path, verbose)

    if filtering:
        filter_propagation(orthofinder_padmet_path, orthofinder_filtered_path, verbose)


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

def addOrthologyInPadmet(orthologue_folder, padmet_folder, output_folder, verbose=False):
    """
    """
    all_orgs = set()
    all_orthologue_files = []
    for _path, _folders, _files in os.walk(orthologue_folder):
        for _folder in _folders:
            all_orgs.add(_folder.replace("Orthologues_","").lower())
        for _file in _files:
            all_orthologue_files.append(os.path.join(_path,_file))
        dict_orthologues = {}
    for org in all_orgs:
        dict_orthologues[org] = dict()
    #dict_orthologue: k = org_id, v = dict: k = gene_id, v = dict: k = org_id, v = set of gene orthologue id
    if verbose:
        print("Extracting orthologues...")
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
                    try:
                        dict_orthologues[org_A.lower()][gene_id_A][org_B.lower()] = set(gene_ids_B)
                    except KeyError:
                        dict_orthologues[org_A.lower()][gene_id_A] = {org_B.lower():set(gene_ids_B)}

    if verbose:
        print("Updating padmets...")
    for padmet_file in [i for i in next(os.walk(padmet_folder))[2]]:
        if verbose:
            print("%s..."%padmet_file)
        org_id = os.path.splitext(padmet_file)[0]
        padmet_path = os.path.join(padmet_folder, padmet_file)
        new_padmet_path = os.path.join(output_folder, "%s.padmet"%org_id)
        org_id = org_id.lower()
        padmet = PadmetSpec(padmet_path)
        for linked_rlt in [rlt for rlt in padmet.getAllRelation() if rlt.type == "is_linked_to"]:
            gene_id = linked_rlt.id_out
            for index, src in enumerate(linked_rlt.misc["SOURCE:ASSIGNMENT"]):
                if 'GENOME:' not in src:
                    ortho_org_id = src.replace("OUTPUT_ORTHOFINDER_FROM_","").lower()
                    ortho_genes_id = ";".join(dict_orthologues[org_id][gene_id][ortho_org_id])
                    new_src = "%s:%s"%(src, ortho_genes_id)
                    linked_rlt.misc["SOURCE:ASSIGNMENT"][index] = new_src
        padmet.generateFile(new_padmet_path)

def filter_propagation(padmet_folder, output_folder, verbose=None):
    propagation_to_remove_file = os.path.join(output_folder, "propagation_to_remove.csv")
    reactions_to_remove_file = os.path.join(output_folder, 'reactions_to_remove.csv')
    padmet_output_folder = output_folder

    if verbose:
        print("Extracting all the relations gene-reaction...")
    dict_rxn_orgs_genes, dict_rxn_ec, dict_orgs_genes = extractRGL(padmet_folder)
    if verbose:
        print("Extracting all the gene propagations...")
    dict_rxn_org_gene_propagation = extractPropagtion(dict_rxn_orgs_genes)
    if verbose:
        print("Writing the file propagation_to_remove...")
    dict_rxn_org_gene_propag_to_remove = extractPropagationToRemove(dict_rxn_org_gene_propagation, output=propagation_to_remove_file)
    if verbose:
        print("Cleaning the Padmet files and writing the reactions_to_remove_file file...")
    cleanPadmet(dict_rxn_org_gene_propag_to_remove, dict_rxn_ec, padmet_folder,
                padmet_output_folder, reactions_to_remove_file)


def extractRGL(padmet_folder):
    """
    extract reactions genes relations.
    It reads all Padmet files in padmet_folder, then it creates three
    dictionaries: dict_rnx_orgs_genes, dict_rnx_ec, dict_org_genes.

    dict_rxn_org_gene: reaction & organism & assigned genes & orthologuous genes.
    dict_rxn_ec[rxn_id] = ec, simple dictionary


    return dict {rxn_id:{org_id:{'FROM-PTOOL': bool, gene_id:set of sources (ex ortho_org_id:ortho_genes)}}}
    """

    all_padmets = [i for i in next(os.walk(padmet_folder))[2]]

    dict_rxn_orgs_genes = {}
    dict_rxn_ec = {}
    dict_org_genes = {}
    for padmet_file in all_padmets:
        padmet_path = os.path.join(padmet_folder, padmet_file)
        org_id = os.path.splitext(padmet_file)[0].upper()
        padmet = PadmetSpec(padmet_path)
        dict_org_genes[org_id] = set()
        for node in padmet.dicOfNode.values():
            if node.type == "reaction":
                rxn_id = node.id
                try:
                    dict_rxn_orgs_genes[rxn_id][org_id] = dict()
                    dict_rxn_ec[rxn_id] = node.misc['EC-NUMBER'][0]
                except KeyError:
                    dict_rxn_orgs_genes[rxn_id] = {org_id:dict()}
                #rxn_id = "PHOSGLYPHOS-RXN"
                from_ptool = False
                for is_linked_rlt in [rlt for rlt in padmet.dicOfRelationIn[rxn_id] if rlt.type == "is_linked_to"]:
                    gene_id = is_linked_rlt.id_out
                    dict_org_genes[org_id].add(gene_id)
                    all_sources = {src.replace("OUTPUT_ORTHOFINDER_FROM_","")  for src in is_linked_rlt.misc["SOURCE:ASSIGNMENT"]}
                    if any (src for src in all_sources if src.startswith("GENOME")):
                        from_ptool = True
                    dict_rxn_orgs_genes[rxn_id][org_id][gene_id] = all_sources
                dict_rxn_orgs_genes[rxn_id][org_id]["FROM-PTOOL"] = from_ptool
    return dict_rxn_orgs_genes, dict_rxn_ec, dict_org_genes


def extractPropagtion(dict_rxn_orgs_genes):
    """
    Propagations are extracted from Padmet networks. Then propagations
    are split bewteen two groups. Those who are orthologues to gene-reaction
    associations that come from Pathway Tools, and those whose associations
    do not come from Pathway Tools.
    """
    dict_rxn_org_gene_propagation = {}
    for rxn_id, rxn_data in dict_rxn_orgs_genes.items():
        #rxn_id ="1.2.1.13-RXN"
        #rxn_data = dict_rxn_orgs_genes[rxn_id]
        dict_rxn_org_gene_propagation[rxn_id] = dict()
        for org_id, org_data in rxn_data.items():
            if org_id not in dict_rxn_org_gene_propagation[rxn_id].keys():
                dict_rxn_org_gene_propagation[rxn_id][org_id] = dict()
            #org_id = 'auxenochlorella_protothecoides'
            #org_data = rxn_data[org_id]
            #because of genes without gene link sources
            if len(org_data.keys()) > 1:
                for gene_id, gene_data in org_data.items():
                    #gene_id, gene_data = list(org_data.items())[2]
                    if gene_id != "FROM-PTOOL":
                        gene_id = gene_id
                        if any (src for src in gene_data if src.startswith("GENOME")):
                            is_from_ptool = True
                        else:
                            is_from_ptool = False
                        for src in gene_data:
                            #src = list(gene_data)[-1]
                            if not src.startswith("GENOME:"):
                                ortho_org_id, ortho_genes_ids = src.split(":")
                                ortho_genes_ids = ortho_genes_ids.split(";")
                                if ortho_org_id not in dict_rxn_org_gene_propagation[rxn_id].keys():
                                    dict_rxn_org_gene_propagation[rxn_id][ortho_org_id] = dict()

                                for ortho_gene_id in ortho_genes_ids:
                                    try:
                                        dict_rxn_org_gene_propagation[rxn_id][org_id][gene_id]["propagation_from_ptool"].add((ortho_org_id, ortho_gene_id))
                                    except KeyError:
                                        dict_rxn_org_gene_propagation[rxn_id][org_id][gene_id] = {"propagation_to_ptool":set(),"propagation_to_not_ptool":set(), "propagation_from_ptool": set([(ortho_org_id, ortho_gene_id)])}

                                    if ortho_gene_id not in dict_rxn_org_gene_propagation[rxn_id][ortho_org_id].keys():
                                        dict_rxn_org_gene_propagation[rxn_id][ortho_org_id][ortho_gene_id] = {"propagation_to_ptool":set(),"propagation_to_not_ptool":set(), "propagation_from_ptool": set()}
                                    if is_from_ptool:
                                        dict_rxn_org_gene_propagation[rxn_id][ortho_org_id][ortho_gene_id]["propagation_to_ptool"].add((org_id, gene_id))
                                    else:
                                        dict_rxn_org_gene_propagation[rxn_id][ortho_org_id][ortho_gene_id]["propagation_to_not_ptool"].add((org_id, gene_id))

    return dict_rxn_org_gene_propagation



def tableExcel(dict_rxn_org_gene_propagation, output):
    with open(output, 'w') as f:
        for rxn_id in ["1.2.1.13-RXN", "1TRANSKETO-RXN", "2TRANSKETO-RXN", "F16ALDOLASE-RXN", "F16BDEPHOS-RXN", "PHOSGLYPHOS-RXN", "PHOSPHORIBULOKINASE-RXN", "RIB5PISOM-RXN", "RIBULOSE-BISPHOSPHATE-CARBOXYLASE-RXN", "RIBULP3EPIM-RXN", "SEDOHEPTULOSE-BISPHOSPHATASE-RXN", "TRIOSEPISOMERIZATION-RXN"]:
            f.write(rxn_id+"\n")
            f.write("\t".join(["gene_id", "prop_to_ptool", "prop_not_to_ptool"])+"\n")
            for org_id, org_data in dict_rxn_org_gene_propagation[rxn_id].items():
                for gene_id, gene_data in org_data.items():
                    if gene_data["propagation_to_ptool"] or gene_data["propagation_to_not_ptool"]:
                        line = "\t".join([gene_id, str(len({i[0] for i in gene_data["propagation_to_ptool"]})), str(len({i[0] for i in gene_data["propagation_to_not_ptool"]}))])+"\n"
                        f.write(line)
            f.write("\n\n")


def extractPropagationToRemove(dict_rxn_org_gene_propagation, output):
    """
    It writes the file propagation_to_remove.
    """
    header = ["reaction_id", "org_id", "gene_id"]
    seuil_ptool = 0

    dict_rxn_org_gene_propag_to_remove = dict()
    with open(output, 'w') as csvfile:
        writer = csv.DictWriter(csvfile, header, delimiter="\t")
        writer.writeheader()
        for rxn_id ,rxn_data in dict_rxn_org_gene_propagation.items():
            nb_org_prop = len(rxn_data.keys())-1
            if nb_org_prop:
                seuil_not_ptool = round(max([(2/nb_org_prop)*100, (5*nb_org_prop)/100]),0)
                for org_id, org_data in rxn_data.items():
                    gene_ids_to_remove = set()
                    for gene_id, gene_data in org_data.items():
                        if len({i[0] for i in gene_data["propagation_to_ptool"]}) <= seuil_ptool:
                            if len({i[0] for i in gene_data["propagation_to_not_ptool"]}) >= seuil_not_ptool:
                                gene_ids_to_remove.add(gene_id)
                    if gene_ids_to_remove:
                        genes_ids = ";".join(gene_ids_to_remove)
                        try:
                            dict_rxn_org_gene_propag_to_remove[rxn_id][org_id] = gene_ids_to_remove
                        except KeyError:
                            dict_rxn_org_gene_propag_to_remove[rxn_id] = {org_id:  gene_ids_to_remove}

                        line = {"reaction_id": rxn_id, "org_id": org_id,
                                "gene_id": genes_ids}
                        writer.writerow(line)
    return dict_rxn_org_gene_propag_to_remove


def cleanPadmet(dict_rxn_org_gene_propag_to_remove, dict_rxn_ec, padmet_folder,
                output_folder, reactions_to_remove_file):
    """
    It cleans the Padmet files and it writes the reactions_to_remove_file file.
    """
    all_padmets = [i for i in next(os.walk(padmet_folder))[2]]

    dict_org_rxn_clean = {}
    for padmet_file in all_padmets:
        #padmet_file = "Thalassiosira_pseudonana.padmet"
        padmet_path = os.path.join(padmet_folder, padmet_file)
        org_id = os.path.splitext(padmet_file)[0].upper()
        padmet = PadmetSpec(padmet_path)
        dict_org_rxn_clean[org_id] = dict()
        for rxn_id in [node.id for node in padmet.dicOfNode.values() if node.type == "reaction" and node.id in dict_rxn_org_gene_propag_to_remove.keys()]:
            dict_org_rxn_clean[org_id][rxn_id] = dict()
            #rxn_id = "1.2.1.13-RXN"
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

    for padmet_file in all_padmets:
        #padmet_file = "Thalassiosira_pseudonana.padmet"
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
                    ec = dict_rxn_ec[rxn_id]
                    gene_id = dict_org_rxn_clean[org_id][rxn_id]
                    line = {"org_id":org_id, "reaction_id":rxn_id,
                            "ec-number":ec, "gene_id":";".join(gene_id)}
                    writer.writerow(line)
                    list_count.append(rxn_id)
    return dict_org_rxn_clean
