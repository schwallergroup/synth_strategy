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
    This function detects if the synthetic route follows a linear synthesis strategy
    (no convergent steps where multiple fragments are combined).
    """
    is_linear = True

    def dfs_traverse(node):
        nonlocal is_linear

        if node["type"] == "reaction":
            try:
                rsmi = node["metadata"]["rsmi"]
                reactants_smiles = rsmi.split(">")[0].split(".")

                # If there are multiple reactants, it's a convergent step
                if len(reactants_smiles) > 1:
                    # Check if all but one are small molecules (potential reagents)
                    non_reagent_count = 0
                    for r in reactants_smiles:
                        mol = Chem.MolFromSmiles(r)
                        if mol and mol.GetNumAtoms() > 5:  # Arbitrary threshold for "non-reagent"
                            non_reagent_count += 1

                    if non_reagent_count > 1:
                        print(f"Convergent step detected: {rsmi}")
                        is_linear = False
            except Exception as e:
                print(f"Error in linear synthesis detection: {e}")

        # Process children
        for child in node.get("children", []):
            dfs_traverse(child)

    # Start traversal
    dfs_traverse(route)

    print(f"Linear synthesis strategy: {is_linear}")
    return is_linear
