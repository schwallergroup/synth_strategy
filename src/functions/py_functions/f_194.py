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
    Detects if the synthesis includes a late-stage nitrile hydrolysis to form an amide.
    """
    # Track if we found the transformation
    found_nitrile_hydrolysis = False

    # SMARTS patterns
    nitrile_pattern = Chem.MolFromSmarts("[#6]C#N")  # Nitrile group
    amide_pattern = Chem.MolFromSmarts("[#6]C(=O)[NH2]")  # Primary amide group

    def dfs_traverse(node, depth=0):
        nonlocal found_nitrile_hydrolysis

        if (
            node["type"] == "reaction" and depth <= 1
        ):  # Check only late-stage reactions (low depth)
            if "rsmi" in node.get("metadata", {}):
                rsmi = node["metadata"]["rsmi"]
                reactants = rsmi.split(">")[0].split(".")
                product = rsmi.split(">")[-1]

                # Check for nitrile hydrolysis
                reactant_mols = [Chem.MolFromSmiles(r) for r in reactants]
                product_mol = Chem.MolFromSmiles(product)

                # Look for nitrile in reactants and amide in product
                has_nitrile_reactant = any(
                    mol and mol.HasSubstructMatch(nitrile_pattern)
                    for mol in reactant_mols
                )
                has_amide_product = product_mol and product_mol.HasSubstructMatch(
                    amide_pattern
                )

                if has_nitrile_reactant and has_amide_product:
                    found_nitrile_hydrolysis = True
                    print(f"Found nitrile hydrolysis at depth {depth}")

        # Traverse children
        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    # Start traversal
    dfs_traverse(route)

    return found_nitrile_hydrolysis
