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
    This function detects if the synthesis follows a linear strategy rather than
    a convergent one, by checking if each reaction involves only one complex fragment.
    """
    is_linear = True

    def dfs_traverse(node):
        nonlocal is_linear

        if node["type"] == "reaction":
            if "rsmi" in node.get("metadata", {}):
                rsmi = node["metadata"]["rsmi"]
                reactants = rsmi.split(">")[0].split(".")

                # Count complex fragments (more than 10 atoms)
                complex_fragments = 0
                for reactant in reactants:
                    mol = Chem.MolFromSmiles(reactant)
                    if mol and mol.GetNumAtoms() > 10:
                        complex_fragments += 1

                if complex_fragments > 1:
                    is_linear = False
                    print(f"Found convergent step with {complex_fragments} complex fragments")

        for child in node.get("children", []):
            dfs_traverse(child)

    dfs_traverse(route)

    if is_linear:
        print("Confirmed linear synthesis strategy")
    return is_linear
