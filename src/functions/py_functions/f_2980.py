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
    Detects if a stereocenter is preserved throughout the synthesis.
    """
    has_stereocenter = False
    stereocenter_preserved = True

    def dfs_traverse(node, depth=0):
        nonlocal has_stereocenter, stereocenter_preserved

        if node["type"] == "mol":
            mol = Chem.MolFromSmiles(node["smiles"])
            if mol:
                # Check if molecule has stereocenters
                chiral_centers = Chem.FindMolChiralCenters(mol, includeUnassigned=False)
                if chiral_centers:
                    has_stereocenter = True
                    print(f"Found stereocenter in molecule at depth {depth}")

        elif node["type"] == "reaction":
            if "rsmi" in node.get("metadata", {}):
                rsmi = node["metadata"]["rsmi"]
                product_smiles = rsmi.split(">")[-1]

                product_mol = (
                    Chem.MolFromSmiles(product_smiles) if product_smiles else None
                )

                if product_mol:
                    # Check if product has stereocenters
                    chiral_centers = Chem.FindMolChiralCenters(
                        product_mol, includeUnassigned=False
                    )
                    if not chiral_centers and has_stereocenter:
                        # If we previously had a stereocenter but now don't, it wasn't preserved
                        stereocenter_preserved = False
                        print(f"Stereocenter not preserved at depth {depth}")

        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    dfs_traverse(route)

    # Return True only if we found stereocenters and they were preserved
    result = has_stereocenter and stereocenter_preserved

    if result:
        print("Stereocenter was preserved throughout the synthesis")

    return result
