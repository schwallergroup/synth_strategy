#!/bin/python

"""LM-defined function for strategy description."""

import copy
import re
from collections import deque

import rdkit
import rdkit.Chem as Chem
from rdkit import Chem
from rdkit.Chem import (
    AllChem,
    Descriptors,
    Lipinski,
    rdChemReactions,
    rdFMCS,
    rdMolDescriptors,
    rdmolops,
)
from rdkit.Chem.Scaffolds import MurckoScaffold


def main(route):
    """
    This function detects the use of TBDMS protection of alcohols in the synthesis.
    """
    tbdms_protection_found = False

    def dfs_traverse(node):
        nonlocal tbdms_protection_found

        if node["type"] == "reaction" and "metadata" in node and "rsmi" in node["metadata"]:
            rsmi = node["metadata"]["rsmi"]
            reactants = rsmi.split(">")[0].split(".")
            product = rsmi.split(">")[-1]

            # Check for TBDMS protection: alcohol + TBDMSCl -> TBDMS ether
            product_mol = Chem.MolFromSmiles(product)
            if product_mol:
                # SMARTS for TBDMS ether: O-Si with tert-butyl and two methyl groups
                tbdms_pattern = Chem.MolFromSmarts("[O][Si]([C])([C])[C]([C])([C])[C]")
                if product_mol.HasSubstructMatch(tbdms_pattern):
                    # Check if any reactant has OH group
                    for reactant in reactants:
                        reactant_mol = Chem.MolFromSmiles(reactant)
                        if reactant_mol:
                            oh_pattern = Chem.MolFromSmarts("[OH]")
                            if reactant_mol.HasSubstructMatch(oh_pattern):
                                tbdms_protection_found = True
                                print("TBDMS protection detected")
                                break

        # Traverse children
        for child in node.get("children", []):
            dfs_traverse(child)

    # Start traversal from the root
    dfs_traverse(route)

    return tbdms_protection_found
