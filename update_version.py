#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
usage:
    update_version.py -f FOLDER

options:
    -h --help     Show help.
    -f <folder>    Path to aucome.
"""

import docopt
import filecmp
import os
import re
import subprocess

def diff_file(file1, file2):
    with open(file1, 'r') as first_file:
        with open(file2, 'r') as second_file:
            diff_f1_f2 = set(first_file)- set(second_file)
    with open(file1, 'r') as first_file:
        with open(file2, 'r') as second_file:
            diff_f2_f1 = set(second_file)- set(first_file)

    return diff_f1_f2, diff_f2_f1

def main():
    args = docopt.docopt(__doc__)

    aucome_folder = args['-f']

    reg_version = r'^\#+VERSION:([0-9.]*)#+'
    first_line = open(aucome_folder+'/'+'release.txt','r').read().split('\n')[0]
    new_version_number = re.match(reg_version,first_line).group(1)

    folder_pathname_new_file = aucome_folder+'/'+'new_file/'
    if not os.path.isdir("{0}".format(folder_pathname_new_file)):
        os.mkdir("{0}".format(folder_pathname_new_file))

    with open(folder_pathname_new_file + 'aucome.py','w') as new_python:
        with open(aucome_folder+'/'+'aucome.py', 'r') as python_script:
            for line in python_script:
                if '__version__ = ' in line:
                    line = line.split(' = ')[0] + ' = "' + str(new_version_number) + '"\n'
                new_python.write(line)

    diff_aucome, diff_aucome_new = diff_file(aucome_folder+'/'+'aucome.py', folder_pathname_new_file + 'aucome.py')

    print('aucome.py')
    print('Old file: ', 'no diff' if diff_aucome == set() else diff_aucome)
    print('New file: ', 'no diff' if diff_aucome_new == set() else diff_aucome_new)
    print('\n')

    with open(folder_pathname_new_file + 'README.md','w') as new_readme:
        with open(aucome_folder+'/'+'README.md', 'r') as old_readme:
            for line in old_readme:
                if 'version: ' in line:
                    line = line.split(': ')[0] + ': ' + str(new_version_number) + '\n'
                new_readme.write(line)

    diff_readme, diff_readme_new = diff_file(aucome_folder+'/'+'README.md', folder_pathname_new_file + 'README.md')

    print('Readme')
    print('Old file: ', 'no diff' if diff_readme == set() else diff_readme)
    print('New file: ', 'no diff' if diff_readme_new == set() else diff_readme_new)
    print('\n')

    with open(folder_pathname_new_file + 'Dockerfile','w') as new_dockerfile:
        with open(aucome_folder+'/'+'Dockerfile', 'r') as old_dockerfile:
            for line in old_dockerfile:
                if 'LABEL Version=' in line:
                    line = line.split('=')[0] + '="' + str(new_version_number) + '"\n'
                new_dockerfile.write(line)

    diff_dockerfile, diff_dockerfile_new = diff_file(aucome_folder+'/'+'Dockerfile', folder_pathname_new_file + 'Dockerfile')

    print('Dockerfile')
    print('Old file: ', 'no diff' if diff_dockerfile == set() else diff_dockerfile)
    print('New file: ', 'no diff' if diff_dockerfile_new == set() else diff_dockerfile_new)

if __name__ == "__main__":
    main()
