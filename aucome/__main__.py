"""
usage:
    aucome --init=ID [-v]
    aucome <command> [<args>...]
    aucome -R --run=ID
    aucome --installPWT=PWT_path [--ptools=ptools_path]
    aucome --uninstallPWT

options:
    -h --help     Show help.
    -R     Give access from container.
    --init=ID    Create folder structure.
    --run=ID    Pathname to the comparison workspace.
    -v     Verbose.

The subcommands are:
    check    Check inputs validity
    reconstruction    Run Pathway Tools
    orthology    Run Orthofinder for crossing orthology between species
    structural    Run check of the structural annotation
    merge    Merge all networks (from Pathway Tools, Orthology and Structural)

    workflow    Run Check, Pathway Tools, Orthofinder and Merging of all networks
    analysis    Analyze results
    compare    Compare group of species

See 'aucome <command> -h' for more information on a specific command.
"""
import configparser
import csv
import docopt
import logging
import mpwt
import os
import pkg_resources
import re
import requests
import shutil
import subprocess
import sys

import aucome

logging.basicConfig(format='%(message)s', level=logging.CRITICAL)
logger = logging.getLogger(__name__)
logging.getLogger("aucome").setLevel(logging.DEBUG)


def main(args=None):
    args = docopt.docopt(__doc__, options_first=True)

    command = args.pop('<command>')
    command_args = args.pop('<args>')


    if args["--init"]:
        run_id = args["--init"]
        create_run(run_id)
        chmod_cmds = ["chmod", "-R", "777", run_id]
        subprocess.call(chmod_cmds)
        return

    # Add permission to all folder in run_id, usefull because all cmd exec from container are root based.
    if args['-R']:
        run_id = args["--run"]
        chmod_cmds = ["chmod", "-R", "777", run_id]
        subprocess.call(chmod_cmds)
        return

    if args['--installPWT']:
        installing_pwt(args['--installPWT'], args['--ptools'])
        return

    if args['--uninstallPWT']:
        uninstalling_pwt()
        return

    if command:
        if command not in ['workflow', 'check', 'reconstruction', 'orthology', 'merge', 'structural', 'analysis', 'compare']:
            sys.exit(command + ' not a valid command: workflow, check, reconstruction, orthology, merge, structural, analysis, compare.')

        if '-h' in command_args:
            getattr(aucome, command).command_help()
            sys.exit()

        # Add command to command_args to be parse by docopt.
        command_args.insert(0,command)

        if command == 'workflow':
            aucome.workflow.workflow_parse_args(command_args)

        elif command == 'check':
            aucome.check.check_parse_args(command_args)

        elif command == 'reconstruction':
            aucome.reconstruction.reconstruction_parse_args(command_args)
        
        elif command == 'orthology':
            aucome.orthology.orthology_parse_args(command_args)

        elif command == 'structural':
            aucome.structural.structural_parse_args(command_args)

        elif command == 'merge':
            aucome.merge.merge_parse_args(command_args)

        elif command == 'analysis':
            aucome.analysis.analysis_parse_args(command_args)

        elif command == 'compare':
            aucome.compare.compare_parse_args(command_args)


def create_run(run_id):
    """
    Create a run folder
    """
    if os.path.isdir('{0}'.format(run_id)):
        print("Run '%s' already exist, remove this folder manually before" %run_id)
    else:
        print('creating Run %s' %run_id)
        os.mkdir('{0}'.format(run_id))
        all_folders = [['studied_organisms'], ['networks'], ['orthology_based'],\
                       ['orthology_based', '0_Orthofinder_WD'], \
                        ['orthology_based', '1_sbml_orthology'], \
                       ['orthology_based', '2_padmet_orthology'], ['orthology_based', '3_padmet_filtered'], ['annotation_based'],\
                       ['annotation_based', 'PGDBs'], ['annotation_based', 'PADMETs'],\
                       ['annotation_based', 'SBMLs'], ['analysis'], ['logs'],\
                       ['networks', 'PADMETs'], ['networks', 'SBMLs'],
                       ['structural_check'], ['structural_check', '0_specifics_reactions'],
                       ['structural_check', '1_blast_results'], ['structural_check', '1_blast_results', 'analysis'], ['structural_check', '1_blast_results', 'tmp'],
                       ['structural_check', '1_blast_results', 'reactions_sequences'],
                       ['structural_check', '2_reactions_to_add'], ['structural_check', '3_PADMETs']]
        for folder in all_folders:
            folder_path = os.path.join(run_id, *folder)
            print('creating folder {0}'.format(folder_path))
            os.mkdir(folder_path)

        group_template_path = os.path.join(run_id, *['analysis', 'group_template.tsv'])
        studied_organisms_path = os.path.join(run_id, 'studied_organisms')
        with open(group_template_path, 'w') as group_file:
            group_writer = csv.writer(group_file, delimiter='\t')
            group_writer.writerow(['all', *os.listdir(studied_organisms_path)])
        config_file_path = os.path.join(run_id, 'config.txt')
        create_config_file(config_file_path, run_id)


def create_config_file(config_file_path, run_id):
    config = configparser.RawConfigParser()
    config.add_section('DATABASE_PATHS')
    config.set('DATABASE_PATHS', 'database_ref_path', '/home/database/BIOCYC/METACYC/23.5/metacyc_23.5.padmet')
    config.add_section('PATHS_IN_RUN')
    config.set('PATHS_IN_RUN', 'run_id', run_id)
    config.set('PATHS_IN_RUN', 'studied_organisms_path', '/studied_organisms')
    config.set('PATHS_IN_RUN', 'orthology_based_path', '/orthology_based')
    config.set('PATHS_IN_RUN', 'orthofinder_wd_path', '/orthology_based/0_Orthofinder_WD')
    config.set('PATHS_IN_RUN', 'orthofinder_sbml_path', '/orthology_based/1_sbml_orthology')
    config.set('PATHS_IN_RUN', 'orthofinder_padmet_path', '/orthology_based/2_padmet_orthology')
    config.set('PATHS_IN_RUN', 'orthofinder_filtered_path', '/orthology_based/3_padmet_filtered')
    config.set('PATHS_IN_RUN', 'annotation_based_path', '/annotation_based')
    config.set('PATHS_IN_RUN', 'pgdb_from_annotation_path', '%(annotation_based_path)s/PGDBs')
    config.set('PATHS_IN_RUN', 'padmet_from_annotation_path', '%(annotation_based_path)s/PADMETs')
    config.set('PATHS_IN_RUN', 'sbml_from_annotation_path', '%(annotation_based_path)s/SBMLs')
    config.set('PATHS_IN_RUN', 'networks_path', '/networks')
    config.set('PATHS_IN_RUN', 'padmet_from_networks_path', '%(networks_path)s/PADMETs')
    config.set('PATHS_IN_RUN', 'sbml_from_networks_path', '%(networks_path)s/SBMLs')
    config.set('PATHS_IN_RUN', 'log_path', '/logs')
    config.set('PATHS_IN_RUN', 'analysis_path', '/analysis')
    config.set('PATHS_IN_RUN', 'analysis_group_file_path', '%(analysis_path)s/group_template.tsv')

    config.set('PATHS_IN_RUN', 'structural_path', '/structural_check')
    config.set('PATHS_IN_RUN', 'structural_specifics_reactions_path', '%(structural_path)s/0_specifics_reactions')
    config.set('PATHS_IN_RUN', 'structural_blast_results_path', '%(structural_path)s/1_blast_results')
    config.set('PATHS_IN_RUN', 'structural_reactions_to_add_path', '%(structural_path)s/2_reactions_to_add')
    config.set('PATHS_IN_RUN', 'structural_padmets_path', '%(structural_path)s/3_PADMETs')
    config.set('PATHS_IN_RUN', 'structural_blast_results_analysis_path', '%(structural_blast_results_path)s/analysis')
    config.set('PATHS_IN_RUN', 'structural_blast_results_tmp_path', '%(structural_blast_results_path)s/tmp')
    config.set('PATHS_IN_RUN', 'structural_blast_results_reactions_sequences_path', '%(structural_blast_results_path)s/reactions_sequences')

    config.add_section('TOOL_PATHS')
    config.set('TOOL_PATHS', 'orthofinder_bin_path', '/programs/OrthoFinder/orthofinder')
    config.add_section('VAR')
    config.set('VAR', 'study_from_annot_prefix', 'output_pathwaytools_')

    # Writing our configuration file to 'example.cfg'
    with open(config_file_path, 'w') as configfile:
        config.write(configfile)


def installing_pwt(pwt_path, input_ptools_local_path):
    """
    Install silently Pathway-Tools in /programs.
    After running this function you need to source the bashrc.
    Args:
        pwt_path (str): Path to the Pathway Tools installer
        input_ptools_local_path (str): Path to local data of Pathway Tools
    """
    ptools_local_path = '/root'
    if input_ptools_local_path:
        if os.path.isdir(input_ptools_local_path):
            ptools_local_path = input_ptools_local_path
        else:
            print(input_ptools_local_path + ' path does not exist, --ptools must be an existing path.')
            return

    cmd_chmods = ['chmod', 'u+x', pwt_path]
    cmd_installs = [pwt_path, '--InstallDir', '/programs/pathway-tools', '--PTOOLS_LOCAL_PATH', ptools_local_path,
                    '--InstallDesktopShortcuts', '0', '--mode', 'unattended']

    print(' '.join(cmd_chmods))
    subprocess.call(cmd_chmods)
    print(' '.join(cmd_installs))
    subprocess.call(cmd_installs)

    # Add Pathway-Tools in the PATH.
    cmd_echo = '''echo 'export PATH="$PATH:/programs/pathway-tools:"' >> ~/.bashrc'''
    print(cmd_echo)
    subprocess.call(cmd_echo, shell=True)

    print('Now you need to source your bash, run:')
    print('source ~/.bashrc')
    return


def uninstalling_pwt():
    """
    Uninstall Pathway-Tools and can delete ptools-local folder.
    """
    def ask_delete_ptools(ptools_path):
        yes_or_no = input('Delete ptools-local folder (y/n)?')
        if yes_or_no == 'y':
            shutil.rmtree(ptools_path)
            print('Uninstallation of Pahtway-Tools and ptools-local done!')
            return
        elif yes_or_no == 'n':
            print('Uninstallation of Pathway-Tools done!.')
            return
        else:
            print('Wrong command')
            ask_delete_ptools(ptools_path)

    ptools_path = mpwt.find_ptools_path()

    cmd_uninstall = ['/programs/pathway-tools/uninstall', '--mode', 'unattended']
    cmd_clean_bash = '''grep -v 'export PATH="$PATH:/programs/pathway-tools:"' ~/.bashrc > ~/temp.bashrc; mv ~/temp.bashrc ~/.bashrc'''

    print(' '.join(cmd_uninstall))
    subprocess.call(cmd_uninstall)

    print(cmd_clean_bash)
    subprocess.call(cmd_clean_bash, shell=True)

    if os.path.isdir('/root/AIC-prefs'):
        shutil.rmtree('/root/AIC-prefs')

    ask_delete_ptools(ptools_path)
    return


if __name__ == "__main__":
    main()
