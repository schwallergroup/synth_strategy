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
    Detects if the synthesis follows a linear strategy without convergent steps
    involving multiple complex fragments.
    """
    is_linear = True

    def dfs_traverse(node, depth=0):
        nonlocal is_linear

        if node["type"] == "reaction":
            if "rsmi" in node.get("metadata", {}):
                rsmi = node["metadata"]["rsmi"]
                reactants_part = rsmi.split(">")[0]

                # Count complex fragments (fragments with more than 10 atoms)
                complex_fragment_count = 0

                for reactant in reactants_part.split("."):
                    try:
                        mol = Chem.MolFromSmiles(reactant)
                        if mol and mol.GetNumAtoms() > 10:
                            complex_fragment_count += 1
                    except:
                        pass

                # If more than one complex fragment, it's a convergent step
                if complex_fragment_count > 1:
                    print(
                        f"Found convergent step at depth {depth} with {complex_fragment_count} complex fragments"
                    )
                    is_linear = False

        # Traverse children
        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    # Start traversal from the root
    dfs_traverse(route)
    return is_linear
