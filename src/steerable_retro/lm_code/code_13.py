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
    This function detects a convergent synthesis strategy where two complex fragments are joined late in the synthesis.
    """
    found_pattern = False

    def dfs_traverse(node, depth=0):
        nonlocal found_pattern

        if node["type"] == "reaction" and depth <= 1:  # Late stage (low depth)
            if "rsmi" in node.get("metadata", {}):
                rsmi = node["metadata"]["rsmi"]
                reactants = rsmi.split(">")[0].split(".")
                product = rsmi.split(">")[-1]

                # Check for complex fragment joining
                reactant_mols = [Chem.MolFromSmiles(r) for r in reactants if Chem.MolFromSmiles(r)]
                product_mol = Chem.MolFromSmiles(product)

                if product_mol and len(reactant_mols) >= 2:
                    # Count complex fragments (>15 heavy atoms)
                    complex_fragments = sum(
                        1 for r in reactant_mols if r and Descriptors.HeavyAtomCount(r) > 15
                    )

                    # Check if product is significantly more complex than individual reactants
                    product_complexity = Descriptors.HeavyAtomCount(product_mol)
                    max_reactant_complexity = max(
                        [Descriptors.HeavyAtomCount(r) for r in reactant_mols if r], default=0
                    )

                    if (
                        complex_fragments >= 2
                        and product_complexity > max_reactant_complexity * 1.3
                    ):
                        print(
                            "Found convergent synthesis with late-stage joining of complex fragments"
                        )
                        found_pattern = True

        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    dfs_traverse(route)
    return found_pattern
