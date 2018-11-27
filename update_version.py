#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
usage:
    update_version.py -f FOLDER

options:
    -h --help     Show help.
    -f <folder>    Path to compare.
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

    compare_folder = args['-f']

    reg_version = r'^\#+VERSION:([0-9.]*)#+'
    first_line = open(compare_folder+'/'+'release.txt','r').read().split('\n')[0]
    new_version_number = re.match(reg_version,first_line).group(1)

    folder_pathname_new_file = compare_folder+'/'+'new_file/'
    if not os.path.isdir("{0}".format(folder_pathname_new_file)):
        os.mkdir("{0}".format(folder_pathname_new_file))

    with open(folder_pathname_new_file + 'compare.py','w') as new_python:
        with open(compare_folder+'/'+'compare.py', 'r') as python_script:
            for line in python_script:
                if '__version__ = ' in line:
                    line = line.split(' = ')[0] + ' = "' + str(new_version_number) + '"\n'
                new_python.write(line)

    diff_compare, diff_compare_new = diff_file(compare_folder+'/'+'compare.py', folder_pathname_new_file + 'compare.py')

    print('Compare')
    print('Old file: ', 'no diff' if diff_compare == set() else diff_compare)
    print('New file: ', 'no diff' if diff_compare_new == set() else diff_compare_new)
    print('\n')

    with open(folder_pathname_new_file + 'README.md','w') as new_readme:
        with open(compare_folder+'/'+'README.md', 'r') as old_readme:
            for line in old_readme:
                if 'version: ' in line:
                    line = line.split(': ')[0] + ': ' + str(new_version_number) + '\n'
                new_readme.write(line)

    diff_readme, diff_readme_new = diff_file(compare_folder+'/'+'README.md', folder_pathname_new_file + 'README.md')

    print('Readme')
    print('Old file: ', 'no diff' if diff_readme == set() else diff_readme)
    print('New file: ', 'no diff' if diff_readme_new == set() else diff_readme_new)
    print('\n')

    with open(folder_pathname_new_file + 'Dockerfile','w') as new_dockerfile:
        with open(compare_folder+'/'+'Dockerfile', 'r') as old_dockerfile:
            for line in old_dockerfile:
                if 'LABEL Version=' in line:
                    line = line.split('=')[0] + '="' + str(new_version_number) + '"\n'
                new_dockerfile.write(line)

    diff_dockerfile, diff_dockerfile_new = diff_file(compare_folder+'/'+'Dockerfile', folder_pathname_new_file + 'Dockerfile')

    print('Dockerfile')
    print('Old file: ', 'no diff' if diff_dockerfile == set() else diff_dockerfile)
    print('New file: ', 'no diff' if diff_dockerfile_new == set() else diff_dockerfile_new)

if __name__ == "__main__":
    main()
