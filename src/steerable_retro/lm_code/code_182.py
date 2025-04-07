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
    This function detects if the synthesis route follows a convergent approach
    where multiple complex fragments are combined.
    """
    max_fragments_combined = 0

    def dfs_traverse(node):
        nonlocal max_fragments_combined

        if node["type"] == "reaction":
            rsmi = node["metadata"]["rsmi"]
            reactants_smiles = rsmi.split(">")[0].split(".")

            # Count non-trivial reactants (more than 5 heavy atoms)
            complex_fragments = 0
            for r in reactants_smiles:
                mol = Chem.MolFromSmiles(r)
                if mol is not None and mol.GetNumHeavyAtoms() > 5:
                    complex_fragments += 1

            if complex_fragments >= 2:
                print(f"Found reaction combining {complex_fragments} complex fragments")
                max_fragments_combined = max(max_fragments_combined, complex_fragments)

        # Traverse children
        for child in node.get("children", []):
            dfs_traverse(child)

    # Start traversal
    dfs_traverse(route)

    is_convergent = max_fragments_combined >= 2
    print(f"Convergent synthesis strategy: {is_convergent}")
    return is_convergent
