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
    This function detects the use of acetylation as a protection strategy.
    It looks for OH → OAc transformations.
    """
    acetylation_count = 0

    def dfs_traverse(node):
        nonlocal acetylation_count

        if node["type"] == "reaction":
            if "rsmi" in node.get("metadata", {}):
                rsmi = node["metadata"]["rsmi"]
                reactants_str = rsmi.split(">")[0]
                product_str = rsmi.split(">")[-1]

                try:
                    reactants_mol = Chem.MolFromSmiles(reactants_str)
                    product_mol = Chem.MolFromSmiles(product_str)

                    # Pattern for alcohol and acetate
                    alcohol_pattern = Chem.MolFromSmarts("[#6]-[#8H1]")
                    acetate_pattern = Chem.MolFromSmarts("[#6]-[#8]-[#6](=[#8])-[#6]")

                    # Check for alcohol to acetate transformation
                    if (
                        reactants_mol
                        and product_mol
                        and reactants_mol.HasSubstructMatch(alcohol_pattern)
                        and product_mol.HasSubstructMatch(acetate_pattern)
                    ):
                        print("Acetylation protection detected")
                        acetylation_count += 1
                except:
                    pass

        # Continue traversal
        for child in node.get("children", []):
            dfs_traverse(child)

    # Start traversal from the root
    dfs_traverse(route)

    # Return True if at least one acetylation is found
    return acetylation_count > 0
