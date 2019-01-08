#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import os
import subprocess

from padmet.classes import PadmetSpec

def test_compare():
    # Set working folder.
    compare_set_working_folder = ['compare', '--setWorkingFolder', os.getcwd()]
    subprocess.call(compare_set_working_folder)
    # Check genbank and create fasta.
    compare_cmd_check = ['compare', '--run=test', '-c']
    subprocess.call(compare_cmd_check)
    # Run Pathway-Tools.
    compare_cmd_pwt = ['compare', '--run=test', '-p']
    subprocess.call(compare_cmd_pwt)
    # Create padmet from PGDB.
    compare_cmd_check = ['compare', '--run=test', '-c']
    subprocess.call(compare_cmd_check)
    # Run Orthofinder.
    compare_cmd_ortho = ['compare', '--run=test', '-o']
    subprocess.call(compare_cmd_ortho)
    # Merge all networks.
    compare_cmd_merge = ['compare', '--run=test', '-d']
    subprocess.call(compare_cmd_merge)

    padmetSpec = PadmetSpec('test/networks/fatty_acid_beta_oxydation_I_1.padmet')
    expected_fabo_reactions = [node.id for node in padmetSpec.dicOfNode.values() if node.type == "reaction"]
    fabo_reactions = ["ACYLCOASYN-RXN", "ACYLCOADEHYDROG-RXN", "ENOYL-COA-DELTA-ISOM-RXN", "ENOYL-COA-HYDRAT-RXN",
                    "OHBUTYRYL-COA-EPIM-RXN", "OHACYL-COA-DEHYDROG-RXN", "KETOACYLCOATHIOL-RXN"]

    assert set(fabo_reactions).issubset(set(expected_fabo_reactions))
