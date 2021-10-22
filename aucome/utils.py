#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Some utils function.
"""

import configparser
import os
import sys


def parse_config_file(run_id):
    config_file_path = "{0}/config.txt".format(run_id)

    if not os.path.exists(config_file_path):
        sys.exit('Error no existing config file: ' + config_file_path)

    config = configparser.ConfigParser()
    config.read(config_file_path)

    #DATABASE_PATHS
    database_path = config.get('DATABASE_PATHS','database_ref_path')
    #PATHS_IN_RUN
    studied_organisms_path = "{0}/{1}".format(run_id, config.get('PATHS_IN_RUN','studied_organisms_path'))
    orthology_based_path = "{0}/{1}".format(run_id, config.get('PATHS_IN_RUN','orthology_based_path'))
    orthofinder_wd_path = "{0}/{1}".format(run_id, config.get('PATHS_IN_RUN','orthofinder_wd_path'))
    orthofinder_sbml_path = "{0}/{1}".format(run_id, config.get('PATHS_IN_RUN','orthofinder_sbml_path'))
    orthofinder_padmet_path = "{0}/{1}".format(run_id, config.get('PATHS_IN_RUN','orthofinder_padmet_path'))
    orthofinder_filtered_path = "{0}/{1}".format(run_id, config.get('PATHS_IN_RUN','orthofinder_filtered_path'))
    annotation_based_path = "{0}/{1}".format(run_id, config.get('PATHS_IN_RUN','annotation_based_path'))
    pgdb_from_annotation_path = "{0}/{1}".format(run_id, config.get('PATHS_IN_RUN','pgdb_from_annotation_path'))
    padmet_from_annotation_path = "{0}/{1}".format(run_id, config.get('PATHS_IN_RUN','padmet_from_annotation_path'))
    sbml_from_annotation_path = "{0}/{1}".format(run_id, config.get('PATHS_IN_RUN','sbml_from_annotation_path'))
    log_path = "{0}/{1}".format(run_id, config.get('PATHS_IN_RUN','log_path'))
    analysis_path = "{0}/{1}".format(run_id, config.get('PATHS_IN_RUN','analysis_path'))
    analysis_group_file_path = "{0}/{1}".format(run_id, config.get('PATHS_IN_RUN','analysis_group_file_path'))
    networks_path = "{0}/{1}".format(run_id, config.get('PATHS_IN_RUN','networks_path'))
    padmet_from_networks_path = "{0}/{1}".format(run_id, config.get('PATHS_IN_RUN','padmet_from_networks_path'))
    sbml_from_networks_path = "{0}/{1}".format(run_id, config.get('PATHS_IN_RUN','sbml_from_networks_path'))
    structural_path = "{0}/{1}".format(run_id, config.get('PATHS_IN_RUN','structural_path'))
    structural_padmets_path = "{0}/{1}".format(run_id, config.get('PATHS_IN_RUN','structural_padmets_path'))

    #TOOL_PATHS
    orthofinder_bin_path = config.get('TOOL_PATHS','orthofinder_bin_path')
    #VAR
    study_from_annot_prefix = config.get('VAR','study_from_annot_prefix')

    config_data = {'database_path': database_path,
                    'studied_organisms_path': studied_organisms_path,
                    'orthology_based_path': orthology_based_path, 'orthofinder_wd_path': orthofinder_wd_path,
                    'orthofinder_sbml_path': orthofinder_sbml_path, 'orthofinder_padmet_path': orthofinder_padmet_path,
                    'orthofinder_filtered_path': orthofinder_filtered_path,
                    'annotation_based_path': annotation_based_path, 'pgdb_from_annotation_path': pgdb_from_annotation_path,
                    'padmet_from_annotation_path': padmet_from_annotation_path, 'sbml_from_annotation_path': sbml_from_annotation_path,
                    'sbml_from_annotation_path': sbml_from_annotation_path, 'log_path': log_path,
                    'networks_path': networks_path, 'orthofinder_bin_path': orthofinder_bin_path,
                    'structural_path': structural_path, 'structural_padmets_path': structural_padmets_path,
                    'study_from_annot_prefix': study_from_annot_prefix, 'analysis_path': analysis_path,
                    'padmet_from_networks_path': padmet_from_networks_path, 'sbml_from_networks_path': sbml_from_networks_path,
                    'analysis_group_file_path': analysis_group_file_path}

    return config_data
