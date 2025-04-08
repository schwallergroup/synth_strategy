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
    This function detects if a synthetic route uses Grignard reagents for C-C bond formation.
    """
    has_grignard = False

    def dfs_traverse(node):
        nonlocal has_grignard

        if node["type"] == "reaction" and "metadata" in node and "rsmi" in node["metadata"]:
            rsmi = node["metadata"]["rsmi"]
            reactants = rsmi.split(">")[0]
            product = rsmi.split(">")[-1]

            # Check for Grignard pattern
            if "MgBr" in reactants or "MgCl" in reactants or "[Mg]" in reactants:
                # Also check for C=O in reactants (typical Grignard target)
                reactant_mol = Chem.MolFromSmiles(reactants)
                if reactant_mol and reactant_mol.HasSubstructMatch(Chem.MolFromSmarts("[#6]=[#8]")):
                    print(f"Found Grignard reaction: {rsmi}")
                    has_grignard = True

        # Traverse children
        for child in node.get("children", []):
            dfs_traverse(child)

    # Start traversal from root
    dfs_traverse(route)

    print(f"Grignard C-C bond formation strategy detected: {has_grignard}")
    return has_grignard
