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
    This function detects if the synthesis follows a linear strategy
    (each reaction has only one product that's used in the next step).

    In a retrosynthetic context, this means each intermediate molecule
    (except the target) is produced by exactly one reaction and used in
    exactly one subsequent reaction.
    """
    # Track molecule usage across the synthesis
    molecule_usage = {}  # smiles -> count of reactions using it

    # Track if we've seen the target molecule
    target_smiles = route["smiles"]

    def dfs_traverse(node, depth=0):
        if node["type"] == "mol":
            # Skip the target molecule (depth 0)
            if depth > 0:
                mol_smiles = node["smiles"]
                if mol_smiles in molecule_usage:
                    molecule_usage[mol_smiles] += 1
                else:
                    molecule_usage[mol_smiles] = 1

                # If a molecule is used in more than one reaction, it's not linear
                if molecule_usage[mol_smiles] > 1:
                    print(f"Molecule {mol_smiles} is used in multiple reactions (non-linear)")
                    return False

        # Continue traversal
        for child in node.get("children", []):
            if not dfs_traverse(child, depth + 1):
                return False

        return True

    # Start traversal from the root
    is_linear = dfs_traverse(route)

    # Debug output
    if is_linear:
        print("Synthesis follows a linear strategy")
    else:
        print("Synthesis does not follow a linear strategy")

    return is_linear
