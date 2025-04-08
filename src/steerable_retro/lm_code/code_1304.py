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
    This function detects if the route contains a fragment coupling strategy
    where two or more complex fragments are combined.
    """
    found = False

    def dfs_traverse(node):
        nonlocal found

        if node["type"] == "reaction" and "metadata" in node and "rsmi" in node["metadata"]:
            rsmi = node["metadata"]["rsmi"]
            reactants_part = rsmi.split(">")[0]
            product = rsmi.split(">")[-1]

            # Count number of fragments in reactants
            reactant_fragments = reactants_part.split(".")

            # Only consider reactions with multiple reactant fragments
            if len(reactant_fragments) >= 2:
                # Check complexity of fragments (simplified: count atoms)
                complex_fragments = 0
                for frag in reactant_fragments:
                    # Count non-hydrogen atoms as a simple complexity measure
                    atom_count = sum(1 for c in frag if c.isupper())
                    if atom_count >= 5:  # Consider fragments with 5+ atoms as complex
                        complex_fragments += 1

                if complex_fragments >= 2:
                    print(f"Found fragment coupling with {complex_fragments} complex fragments")
                    found = True

        for child in node.get("children", []):
            dfs_traverse(child)

    dfs_traverse(route)

    return found
