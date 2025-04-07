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
    This function detects if the synthesis follows a linear strategy (as opposed to convergent).
    """
    max_fragments_per_reaction = 0

    def dfs_traverse(node):
        nonlocal max_fragments_per_reaction

        if node["type"] == "reaction" and "rsmi" in node.get("metadata", {}):
            rsmi = node["metadata"]["rsmi"]
            reactants_part = rsmi.split(">")[0]

            # Count the number of reactant fragments
            reactants = reactants_part.split(".")
            num_fragments = len(reactants)

            # Count only "complex" fragments (more than 10 atoms)
            complex_fragments = 0
            for reactant in reactants:
                mol = Chem.MolFromSmiles(reactant)
                if mol and mol.GetNumAtoms() > 10:
                    complex_fragments += 1

            max_fragments_per_reaction = max(max_fragments_per_reaction, complex_fragments)

        for child in node.get("children", []):
            dfs_traverse(child)

    dfs_traverse(route)

    # Linear synthesis typically has at most one complex fragment per reaction
    result = max_fragments_per_reaction <= 1

    print(
        f"Linear synthesis strategy: {result} (max complex fragments: {max_fragments_per_reaction})"
    )
    return result
