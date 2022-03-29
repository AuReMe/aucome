#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import os
import subprocess

from padmet.classes import PadmetSpec

def test_aucome():
    # Set working folder.
    aucome_set_working_folder = ['aucome', '--setWorkingFolder', os.getcwd()]
    subprocess.call(aucome_set_working_folder)
    # Check genbank and create fasta.
    aucome_cmd_check = ['aucome', 'check', '--run=test']
    subprocess.call(aucome_cmd_check)
    # Run Pathway-Tools.
    aucome_cmd_pwt = ['aucome', 'reconstruction', '--run=test']
    subprocess.call(aucome_cmd_pwt)
    # Run Orthofinder.
    aucome_cmd_ortho = ['aucome', 'orthology', '--run=test']
    subprocess.call(aucome_cmd_ortho)
    # Spontaneous all networks.
    aucome_cmd_spontaneous = ['aucome', 'draft', '--run=test']
    subprocess.call(aucome_cmd_spontaneous)

    padmetSpec = PadmetSpec('test/networks/fatty_acid_beta_oxydation_I_1.padmet')
    expected_fabo_reactions = [node.id for node in padmetSpec.dicOfNode.values() if node.type == "reaction"]
    fabo_reactions = ["ACYLCOASYN-RXN", "ACYLCOADEHYDROG-RXN", "ENOYL-COA-DELTA-ISOM-RXN", "ENOYL-COA-HYDRAT-RXN",
                    "OHBUTYRYL-COA-EPIM-RXN", "OHACYL-COA-DEHYDROG-RXN", "KETOACYLCOATHIOL-RXN"]

    assert set(fabo_reactions).issubset(set(expected_fabo_reactions))
