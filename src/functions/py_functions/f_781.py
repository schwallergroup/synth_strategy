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
    This function detects if the synthetic route follows a linear synthesis pattern
    with a late-stage coupling of two major fragments.
    """
    late_coupling_found = False

    def dfs_traverse(node, depth=0):
        nonlocal late_coupling_found

        if node["type"] == "reaction" and depth <= 1:  # Focus on late-stage reactions
            rsmi = node["metadata"]["rsmi"]
            reactants_smiles = rsmi.split(">")[0].split(".")

            # Check if we have multiple complex reactants (indicating fragment coupling)
            complex_reactants = 0
            for r_smiles in reactants_smiles:
                r = Chem.MolFromSmiles(r_smiles)
                if (
                    r and r.GetNumAtoms() > 5
                ):  # Consider reactants with >5 atoms as complex
                    complex_reactants += 1

            if complex_reactants >= 2:
                late_coupling_found = True
                print(
                    f"Late-stage coupling of multiple complex fragments detected at depth {depth}"
                )

        # Traverse children
        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    # Start traversal
    dfs_traverse(route)

    print(f"Linear synthesis with late-stage coupling detected: {late_coupling_found}")
    return late_coupling_found
