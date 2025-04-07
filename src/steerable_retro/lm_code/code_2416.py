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
    This function detects if the synthesis follows a linear fragment coupling strategy
    with sequential addition of fragments.
    """
    fragment_couplings = 0
    max_depth = 0

    def dfs_traverse(node, depth=0):
        nonlocal fragment_couplings, max_depth
        max_depth = max(max_depth, depth)

        if node["type"] == "reaction":
            # Extract reactants and product
            rsmi = node["metadata"]["rsmi"]
            reactants_smiles = rsmi.split(">")[0].split(".")

            # If there are multiple reactants, it might be a fragment coupling
            if len(reactants_smiles) >= 2:
                # Check if reactants are substantial fragments (not just small reagents)
                reactant_mols = [Chem.MolFromSmiles(smi) for smi in reactants_smiles]
                substantial_fragments = 0

                for mol in reactant_mols:
                    if (
                        mol and mol.GetNumHeavyAtoms() > 6
                    ):  # Consider fragments with >6 atoms as substantial
                        substantial_fragments += 1

                if substantial_fragments >= 2:
                    print(f"Fragment coupling detected at depth {depth}")
                    fragment_couplings += 1

        # Traverse children
        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    # Start traversal
    dfs_traverse(route)

    # Check if we have multiple fragment couplings and they occur sequentially
    is_linear_strategy = fragment_couplings >= 2 and fragment_couplings <= max_depth

    if is_linear_strategy:
        print(f"Linear fragment coupling strategy detected with {fragment_couplings} couplings")

    return is_linear_strategy
