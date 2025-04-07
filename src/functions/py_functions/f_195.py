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
    Detects if the synthesis includes reduction of a nitro group to an aniline.
    """
    # Track if we found the transformation
    found_nitro_reduction = False

    # SMARTS patterns
    nitro_pattern = Chem.MolFromSmarts("[#6]-[N+](=[O])[O-]")  # Nitro group
    aniline_pattern = Chem.MolFromSmarts("[#6]-[NH2]")  # Aniline group

    def dfs_traverse(node, depth=0):
        nonlocal found_nitro_reduction

        if node["type"] == "reaction":
            if "rsmi" in node.get("metadata", {}):
                rsmi = node["metadata"]["rsmi"]
                reactants = rsmi.split(">")[0].split(".")
                product = rsmi.split(">")[-1]

                # Check for nitro reduction
                reactant_mols = [Chem.MolFromSmiles(r) for r in reactants]
                product_mol = Chem.MolFromSmiles(product)

                # Look for nitro in reactants and aniline in product
                has_nitro_reactant = any(
                    mol and mol.HasSubstructMatch(nitro_pattern)
                    for mol in reactant_mols
                )
                has_aniline_product = product_mol and product_mol.HasSubstructMatch(
                    aniline_pattern
                )

                if has_nitro_reactant and has_aniline_product:
                    found_nitro_reduction = True
                    print(f"Found nitro reduction at depth {depth}")

        # Traverse children
        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    # Start traversal
    dfs_traverse(route)

    return found_nitro_reduction
