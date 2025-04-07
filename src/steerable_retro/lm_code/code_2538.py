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
    This function detects a linear synthesis strategy with sequential amide bond formations.
    """
    # Track amide formations at different depths
    amide_formations = []

    # SMARTS pattern for amide bond
    amide_pattern = Chem.MolFromSmarts("[#6](=[#8])-[#7]")

    def dfs_traverse(node, depth=0):
        if node["type"] == "reaction":
            # Extract reactants and product
            rsmi = node["metadata"]["rsmi"]
            reactants_smiles = rsmi.split(">")[0].split(".")
            product_smiles = rsmi.split(">")[-1]

            # Count amide bonds in reactants and product
            reactants_amide_count = sum(
                len(Chem.MolFromSmiles(r).GetSubstructMatches(amide_pattern))
                for r in reactants_smiles
                if Chem.MolFromSmiles(r) is not None
            )

            product_mol = Chem.MolFromSmiles(product_smiles)
            if product_mol is not None:
                product_amide_count = len(product_mol.GetSubstructMatches(amide_pattern))

                # If product has more amide bonds than reactants, an amide formation occurred
                if product_amide_count > reactants_amide_count:
                    amide_formations.append(depth)
                    print(f"Amide bond formation detected at depth {depth}")

        # Traverse children
        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    # Start traversal
    dfs_traverse(route)

    # Check if we have multiple amide formations at different depths
    if len(amide_formations) >= 2:
        print(f"Linear amide formation strategy detected at depths: {sorted(amide_formations)}")
        return True
    return False
