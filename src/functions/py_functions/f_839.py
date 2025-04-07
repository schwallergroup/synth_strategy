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
    Detects if the synthesis uses a convergent approach with multiple amide couplings
    to assemble fragments.
    """
    amide_couplings = []
    fragments = set()

    def dfs_traverse(node, depth=0):
        nonlocal amide_couplings, fragments

        if node["type"] == "reaction":
            # Extract reactants and product
            rsmi = node["metadata"]["rsmi"]
            reactants_smiles = rsmi.split(">")[0].split(".")
            product_smiles = rsmi.split(">")[-1]

            product_mol = Chem.MolFromSmiles(product_smiles)

            # Check for amide formation
            amide_pattern = Chem.MolFromSmarts("[NX3][CX3](=[OX1])")
            if product_mol.HasSubstructMatch(amide_pattern):
                # Count number of fragments being combined
                if len(reactants_smiles) >= 2:
                    amide_couplings.append((depth, len(reactants_smiles)))
                    for r in reactants_smiles:
                        fragments.add(r)
                    print(
                        f"Found amide coupling at depth {depth} combining {len(reactants_smiles)} fragments"
                    )

        # Traverse children
        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    # Start traversal from root
    dfs_traverse(route)

    # Check if the strategy is present
    strategy_present = len(amide_couplings) >= 2 and len(fragments) >= 3

    print(f"Convergent amide coupling strategy detected: {strategy_present}")
    print(f"Number of amide couplings: {len(amide_couplings)}")
    print(f"Number of unique fragments: {len(fragments)}")

    return strategy_present
