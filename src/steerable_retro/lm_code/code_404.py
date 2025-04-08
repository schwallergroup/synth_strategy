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
    Detects a convergent synthesis strategy where two or more complex fragments
    are joined in a late-stage reaction.
    """
    found_convergent_step = False

    def dfs_traverse(node, depth=0):
        nonlocal found_convergent_step

        if node["type"] == "reaction" and depth <= 1:  # Focus on late-stage reactions
            if "metadata" in node and "rsmi" in node["metadata"]:
                rsmi = node["metadata"]["rsmi"]
                reactants_str = rsmi.split(">")[0]

                # Count number of distinct reactants
                reactants = reactants_str.split(".")

                # Filter out small molecules/reagents (less than 10 atoms)
                complex_reactants = []
                for reactant in reactants:
                    try:
                        mol = Chem.MolFromSmiles(reactant)
                        if mol and mol.GetNumAtoms() > 10:  # Consider only substantial fragments
                            complex_reactants.append(reactant)
                    except:
                        continue

                if len(complex_reactants) >= 2:
                    found_convergent_step = True
                    print(
                        f"Found convergent step at depth {depth} with {len(complex_reactants)} complex fragments"
                    )

        # Traverse children
        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    # Start traversal from the root
    dfs_traverse(route)

    return found_convergent_step
