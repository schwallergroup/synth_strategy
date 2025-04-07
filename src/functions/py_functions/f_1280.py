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
    This function detects a synthetic strategy involving selective mono-functionalization
    of a symmetrical diol.
    """
    has_selective_diol_functionalization = False

    def dfs_traverse(node):
        nonlocal has_selective_diol_functionalization

        if node["type"] == "reaction" and "rsmi" in node.get("metadata", {}):
            rsmi = node["metadata"]["rsmi"]
            reactants = rsmi.split(">")[0].split(".")
            product = rsmi.split(">")[-1]

            # Check for diol in reactants
            diol_pattern = Chem.MolFromSmarts("[OX2H][CX4][CX4][OX2H]")

            # Check for mono-protected product (one OH, one protected OH)
            protected_pattern = Chem.MolFromSmarts("[OX2H][CX4][CX4][OX2][Si]")

            for r in reactants:
                try:
                    mol = Chem.MolFromSmiles(r)
                    if mol and mol.HasSubstructMatch(diol_pattern):
                        prod_mol = Chem.MolFromSmiles(product)
                        if prod_mol and prod_mol.HasSubstructMatch(protected_pattern):
                            has_selective_diol_functionalization = True
                            print("Found selective diol functionalization")
                except:
                    continue

        # Recursively traverse children
        for child in node.get("children", []):
            dfs_traverse(child)

    # Start traversal from the root
    dfs_traverse(route)

    return has_selective_diol_functionalization
