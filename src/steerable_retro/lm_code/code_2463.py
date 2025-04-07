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
    Detects amide formation from acid chloride, a common strategy for introducing amide functionality.
    """
    acid_chloride_to_amide = False

    def dfs_traverse(node, depth=0):
        nonlocal acid_chloride_to_amide

        if node["type"] == "reaction":
            if "rsmi" in node.get("metadata", {}):
                rsmi = node["metadata"]["rsmi"]
                reactants = rsmi.split(">")[0].split(".")
                product = rsmi.split(">")[-1]

                # Check for acid chloride pattern in reactants
                acid_chloride_pattern = Chem.MolFromSmarts("[#6](=[#8])[Cl]")

                # Check for amide pattern in product
                amide_pattern = Chem.MolFromSmarts("[#6](=[#8])[#7]")

                has_acid_chloride = any(
                    Chem.MolFromSmiles(r) is not None
                    and Chem.MolFromSmiles(r).HasSubstructMatch(acid_chloride_pattern)
                    for r in reactants
                )

                product_mol = Chem.MolFromSmiles(product)
                has_amide = product_mol is not None and product_mol.HasSubstructMatch(amide_pattern)

                if has_acid_chloride and has_amide:
                    acid_chloride_to_amide = True
                    print(f"Found acid chloride to amide conversion at depth {depth}")

        # Continue traversal
        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    # Start traversal
    dfs_traverse(route)

    return acid_chloride_to_amide
