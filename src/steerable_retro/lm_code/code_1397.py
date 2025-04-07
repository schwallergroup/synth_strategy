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
    Detects if the synthesis route uses a convergent approach with multiple fragment couplings.
    Counts reactions where 2+ fragments are combined.
    """
    fragment_coupling_count = 0

    def dfs_traverse(node):
        nonlocal fragment_coupling_count

        if node["type"] == "reaction":
            # Extract reactants and product
            rsmi = node["metadata"]["rsmi"]
            reactants_smiles = rsmi.split(">")[0].split(".")

            # If reaction has 2+ reactants, it's potentially a fragment coupling
            if len(reactants_smiles) >= 2:
                # Check if reactants are significant fragments (not just reagents)
                significant_fragments = 0
                for smi in reactants_smiles:
                    mol = Chem.MolFromSmiles(smi)
                    if mol:
                        # Consider a fragment significant if it has at least one ring
                        if Chem.GetSSSR(mol):
                            significant_fragments += 1

                if significant_fragments >= 2:
                    fragment_coupling_count += 1
                    print(f"Fragment coupling detected in reaction: {rsmi}")

        # Traverse children
        for child in node.get("children", []):
            dfs_traverse(child)

    # Start traversal
    dfs_traverse(route)

    # Consider it a convergent strategy if there are at least 2 fragment couplings
    is_convergent = fragment_coupling_count >= 2
    print(
        f"Convergent synthesis strategy: {is_convergent} (fragment couplings: {fragment_coupling_count})"
    )
    return is_convergent
