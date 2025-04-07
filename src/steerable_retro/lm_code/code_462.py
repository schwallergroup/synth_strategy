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
    Detects if the synthetic route includes a cyanation reaction (conversion of aryl halide to nitrile).
    """
    has_cyanation = False

    def dfs_traverse(node):
        nonlocal has_cyanation

        if node["type"] == "reaction":
            if "rsmi" in node["metadata"]:
                rsmi = node["metadata"]["rsmi"]
                reactants_smiles = rsmi.split(">")[0]
                product_smiles = rsmi.split(">")[-1]

                # Check for aryl halide in reactants
                aryl_halide_pattern = Chem.MolFromSmarts("[c]-[#9,#17,#35,#53]")
                reactants_mol = Chem.MolFromSmiles(reactants_smiles)

                # Check for nitrile in product
                nitrile_pattern = Chem.MolFromSmarts("[C]#[N]")
                product_mol = Chem.MolFromSmiles(product_smiles)

                if reactants_mol and product_mol:
                    if reactants_mol.HasSubstructMatch(
                        aryl_halide_pattern
                    ) and product_mol.HasSubstructMatch(nitrile_pattern):
                        has_cyanation = True
                        print(f"Found cyanation reaction: {rsmi}")

        # Traverse children
        for child in node.get("children", []):
            dfs_traverse(child)

    # Start traversal from the root
    dfs_traverse(route)

    return has_cyanation
