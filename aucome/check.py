#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
usage:
    aucome check --run=ID [--cpu=INT] [-v]

options:
    --run=ID    Pathname to the comparison workspace.
    --cpu=INT     Number of cpu to use for the multiprocessing (if none use 1 cpu).
    -v     Verbose.

"""

import csv
import docopt
import os

from padmet.utils.connection import gbk_to_faa, pgdb_to_padmet, sbmlGenerator

from aucome.utils import parse_config_file

from Bio import SeqIO
from multiprocessing import Pool


def command_help():
    print(docopt.docopt(__doc__))


def check_parse_args(command_args):
    args = docopt.docopt(__doc__, argv=command_args)
    run_id = args['--run']
    verbose = args['-v']

    if args["--cpu"]:
        nb_cpu_to_use = int(args["--cpu"])
    else:
        nb_cpu_to_use = 1

    run_check(run_id, nb_cpu_to_use, verbose)


def run_check(run_id, nb_cpu_to_use, verbose):

    config_data = parse_config_file(run_id)

    padmet_from_annotation_path = config_data['padmet_from_annotation_path']
    study_from_annot_prefix = config_data['study_from_annot_prefix']
    sbml_from_annotation_path = config_data['sbml_from_annotation_path']
    database_path = config_data['database_path']
    pgdb_from_annotation_path = config_data['pgdb_from_annotation_path']
    studied_organisms_path = config_data['studied_organisms_path']
    model_organisms_path = config_data['model_organisms_path']
    analysis_group_file_path = config_data['analysis_group_file_path']

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

    # Update group file in analysis
    if not os.path.exists(analysis_group_file_path):
        with open(analysis_group_file_path, 'w') as group_file:
            group_writer = csv.writer(group_file, delimiter='\t')
            group_writer.writerow(['all', *all_study_name])
    else:
        groups_data = []
        with open(analysis_group_file_path, 'r') as group_file:
            group_reader = csv.reader(group_file, delimiter='\t')
            for row in group_reader:
                groups = [org_name for org_name in row[1:] if org_name]
                groups_data.append((row[0], groups))

        # Check if 'all' row matches species in study_organisms.
        if sorted(groups_data[0][1]) != sorted(all_study_name):
            with open(analysis_group_file_path, 'w') as group_file:
                group_writer = csv.writer(group_file, delimiter='\t')
                group_writer.writerow(['all', *all_study_name])
                for group in groups_data:
                    if group[0] != 'all':
                        group_writer.writerow([group[0], *group[1]])

    aucome_pool = Pool(nb_cpu_to_use)

    if verbose:
        print('Checking genbank file.')
    study_faa_data = []
    for study_name in all_study_name:
        faa_path = "{0}/{1}/{1}.faa".format(studied_organisms_path, study_name)
        tmp_faa_data = {'study_name': study_name, 'faa_path': faa_path, 'gbk_file': all_study_gbk[study_name],
                        'studied_organisms_path': studied_organisms_path,
                        'verbose': verbose}
        study_faa_data.append(tmp_faa_data)
    aucome_pool.map(check_create_faa, study_faa_data)

    #k = folder_name in studied_org_path, v = path to faa in this folder, faa name should be folder_name.faa
    all_study_faa = dict([(study_name, "{0}/{1}/{1}.faa".format(studied_organisms_path, study_name))
                          if os.path.isfile("{0}/{1}/{1}.faa".format(studied_organisms_path, study_name))
                          else (study_name, '')
                          for study_name in all_study_name])

    study_model_data = []
    for model_name in all_model_name:
        faa_path = "{0}/{1}/{1}.faa".format(model_organisms_path, model_name)
        tmp_model_data = {'model_name': model_name, 'faa_path': faa_path, 'gbk_file': all_model_gbk[model_name],
                            'verbose': verbose}
        study_model_data.append(tmp_model_data)
    aucome_pool.map(create_faa_model, study_model_data)

    #k = folder_name in model_organisms_path, v = path to faa in this folder, faa name should be folder_name.faa
    all_model_faa = dict([(model_name, "{0}/{1}/{1}.faa".format(model_organisms_path, model_name))
                          if os.path.isfile("{0}/{1}/{1}.faa".format(model_organisms_path, model_name))
                          else (model_name, '')
                          for model_name in all_model_name])

    study_padmet_data = []
    for study_name in all_study_name:
        padmet_file = "{0}/{1}{2}.padmet".format(padmet_from_annotation_path, study_from_annot_prefix, study_name)
        pgdb_folder = all_study_pgdb[study_name]
        tmp_padmet_data = {'study_name': study_name, 'pgdb_folder': pgdb_folder, 'verbose': verbose, 
                           'padmet_file': padmet_file, 'database_path': database_path}
        study_padmet_data.append(tmp_padmet_data)
    aucome_pool.map(create_padmet_from_pgdb, study_padmet_data)

    all_study_padmet = dict([(study_name, "{0}/{1}{2}.padmet".format(padmet_from_annotation_path, study_from_annot_prefix, study_name))
                          if os.path.isfile("{0}/{1}{2}.padmet".format(padmet_from_annotation_path, study_from_annot_prefix, study_name))
                          else (study_name, '')
                          for study_name in all_study_name])

    study_sbml_data = []
    for study_name in all_study_padmet:
        sbml_file = "{0}/{1}{2}.sbml".format(sbml_from_annotation_path, study_from_annot_prefix, study_name)
        padmet_file = all_study_padmet[study_name]
        tmp_sbml_data = {'sbml_file': sbml_file, 'padmet_file': padmet_file,
                         'study_name': study_name, 'verbose': verbose}
        study_sbml_data.append(tmp_sbml_data)
    aucome_pool.map(create_sbml, study_sbml_data)

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

    aucome_pool.close()
    aucome_pool.join()


def check_create_faa(tmp_faa_data):
    study_name = tmp_faa_data['study_name']
    faa_path = tmp_faa_data['faa_path']
    gbk_file = tmp_faa_data['gbk_file']
    verbose = tmp_faa_data['verbose']
    studied_organisms_path = tmp_faa_data['studied_organisms_path']
    checking_genbank(study_name, studied_organisms_path, verbose)

    #create Faa from gbk if no faa found
    if not os.path.isfile(faa_path) and gbk_file:
        if verbose:
            print("Creating faa from gbk for %s" %study_name)
        gbk_to_faa.gbk_to_faa(gbk_file=gbk_file, output=faa_path, verbose=verbose)


def create_faa_model(tmp_model_data):
    model_name = tmp_model_data['model_name']
    gbk_file = tmp_model_data['gbk_file']
    verbose = tmp_model_data['verbose']
    faa_path = tmp_model_data['faa_path']

    if not os.path.isfile(faa_path) and gbk_file:
        if verbose:
            print("Creating faa from gbk for %s" %model_name)
        gbk_to_faa.gbk_to_faa(gbk_file=gbk_file, output=faa_path, verbose=verbose)


def create_padmet_from_pgdb(tmp_padmet_data):
    study_name = tmp_padmet_data['study_name']
    pgdb_folder = tmp_padmet_data['pgdb_folder']
    verbose = tmp_padmet_data['verbose']
    padmet_file = tmp_padmet_data['padmet_file']
    database_path = tmp_padmet_data['database_path']

    if not os.path.isfile(padmet_file) and pgdb_folder:
        if verbose:
            print("Creating padmet from pgdb for %s" %study_name)
        pgdb_to_padmet.from_pgdb_to_padmet(pgdb_folder=pgdb_folder, output=padmet_file, padmetRef_file=database_path, source="genome", extract_gene=True, no_orphan=True, verbose=verbose)



def create_sbml(tmp_sbml_data):
    sbml_file = tmp_sbml_data['sbml_file']
    padmet_file = tmp_sbml_data['padmet_file']
    study_name = tmp_sbml_data['study_name']
    verbose = tmp_sbml_data['verbose']

    if not os.path.isfile(sbml_file) and padmet_file:
        if verbose:
            print("Creating sbml from padmet for %s" %study_name)
        sbmlGenerator.padmet_to_sbml(padmet_file=padmet_file, output=sbml_file, sbml_lvl=3, verbose=True)



def checking_genbank(genbank_file_name, studied_organisms_path, verbose):
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
        fix_genbank_file(genbank_file_name, fix_name, fix_dot_protein_seq, studied_organisms_path, verbose)


def adapt_gene_id(gene_id, longest_gene_number_length):
    """
    Input: a gene ID like g_1 and the longest_gene_number_length (5 if you have a gene like g_12356).
    Return a new gene ID with a completion of 0 like: g_00001.
    """
    gene_prefix = '_'.join(gene_id.split('_')[:-1])
    gene_number = gene_id.split('_')[-1]

    new_gene_id = gene_prefix + '_' + gene_number.zfill(longest_gene_number_length)

    return new_gene_id


def fix_genbank_file(genbank_file_name, fix_name, fix_dot_protein_seq, studied_organisms_path, verbose):
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