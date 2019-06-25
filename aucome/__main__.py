"""
usage:
    aucome --init=ID [-v]
    aucome <command> [<args>...]
    aucome -R --run=ID
    aucome --version
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
    reconstruction    Run Pathway-Tools
    orthology    Run Orthofinder
    draft    Merge all networks
    workflow    Run Orthofinder, Pathway and merge all networks

See 'aucome <command> -h' for more information on a specific command.
"""
import configparser
import docopt
import eventlet
import mpwt
import os
import re
import requests
import subprocess
import sys

import aucome

__version__ = "0.4"
release_on_gitlab = "https://gitlab.inria.fr/DYLISS/compare_metabo/raw/master/release.txt"


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

    if args["--version"]:
        online_version = get_version()
        current_version = __version__
        if online_version:
            print("You are using the version %s, the latest is %s" %(current_version, online_version))
        else:
            print('No internet connection. Skip checking AuCoMe version.')
        return

    if command:
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

        elif command == 'draft':
            aucome.draft.draft_parse_args(command_args)

        if '-h' in args:
            getattr(aucome, command).command_help()


def create_run(run_id):
    """
    Create a run folder
    """
    if os.path.isdir("{0}".format(run_id)):
        print("Run '%s' already exist, remove this folder manually before" %run_id)
    else:
        print("creating Run %s" %run_id)
        os.mkdir("{0}".format(run_id))
        all_folders = ["studied_organisms", "model_organisms", "networks", "orthology_based",\
                       "orthology_based/Orthofinder_WD", "annotation_based",\
                       "annotation_based/PGDBs", "annotation_based/PADMETs",\
                       "annotation_based/SBMLs", "analysis", "logs"]
        for folder in all_folders:
            print("creating folder {0}/{1}".format(run_id, folder))
            os.mkdir("{0}/{1}".format(run_id, folder))
        config_file_path = "{0}/config.txt".format(run_id)
        create_config_file(config_file_path, run_id)


def create_config_file(config_file_path, run_id):
    config = configparser.RawConfigParser()
    config.add_section('DATABASE_PATHS')
    config.set('DATABASE_PATHS', 'database_ref_path', '/home/database/BIOCYC/METACYC/22.0_enhanced/metacyc_22.0_enhanced.padmet')
    config.set('DATABASE_PATHS', 'mnx_rxn_path', '/home/database/MNX/reac_xref.tsv')
    config.set('DATABASE_PATHS', 'mnx_cpd_path', '/home/database/MNX/chem_xref.tsv')
    config.add_section('PATHS_IN_RUN')
    config.set('PATHS_IN_RUN', 'run_id', run_id)
    config.set('PATHS_IN_RUN', 'studied_organisms_path', '/studied_organisms')
    config.set('PATHS_IN_RUN', 'model_organisms_path', '/model_organisms')
    config.set('PATHS_IN_RUN', 'orthology_based_path', '/orthology_based')
    config.set('PATHS_IN_RUN', 'orthofinder_wd_path', '/orthology_based/Orthofinder_WD')
    config.set('PATHS_IN_RUN', 'annotation_based_path', '/annotation_based')
    config.set('PATHS_IN_RUN', 'pgdb_from_annotation_path', '%(annotation_based_path)s/PGDBs')
    config.set('PATHS_IN_RUN', 'padmet_from_annotation_path', '%(annotation_based_path)s/PADMETs')
    config.set('PATHS_IN_RUN', 'sbml_from_annotation_path', '%(annotation_based_path)s/SBMLs')
    config.set('PATHS_IN_RUN', 'networks_path', '/networks')
    config.set('PATHS_IN_RUN', 'log_path', '/logs')
    config.add_section('TOOL_PATHS')
    config.set('TOOL_PATHS', 'orthofinder_bin_path', '/programs/OrthoFinder-2.3.3/orthofinder')
    config.set('TOOL_PATHS', 'padmet_utils_path', '/programs/padmet-utils')
    config.add_section('VAR')
    config.set('VAR', 'study_from_annot_prefix', 'output_pathwaytools_')

    # Writing our configuration file to 'example.cfg'
    with open(config_file_path, 'w') as configfile:
        config.write(configfile)

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

def main_workflow(run_id, orthofidner_sequence_tool, orthogroups_option, cpu_number, verbose, log):
    a = 411
def main_check(run_id, cpu_number, verbose):
    a = 411
def main_reconstruction(run_id, number_cpu, verbose, log_folder):
    a = 411
def main_orthology(run_id, orthofidner_sequence_tool, orthogroups_option, cpu_number, verbose):
    a = 411

if __name__ == "__main__":
    main()