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
    Detects a convergent synthesis strategy where two complex fragments
    (>8 atoms each) are combined in the final step.
    """
    is_convergent = False

    def dfs_traverse(node, depth=0):
        nonlocal is_convergent

        # Check if this is a reaction node
        if node["type"] == "reaction":
            # Get reaction SMILES
            rsmi = node.get("metadata", {}).get("rsmi", "")
            if rsmi:
                try:
                    # Extract reactants
                    reactants_part = rsmi.split(">")[0]
                    reactants = reactants_part.split(".")

                    # Count complex fragments (>8 atoms)
                    complex_fragments = 0

                    for reactant in reactants:
                        if not reactant:
                            continue

                        r_mol = Chem.MolFromSmiles(reactant)
                        if not r_mol:
                            continue

                        if r_mol.GetNumAtoms() > 8:
                            complex_fragments += 1

                    # If this is a final or late-stage reaction (depth <= 1) with multiple complex fragments
                    if depth <= 1 and complex_fragments >= 2:
                        is_convergent = True
                        print(
                            f"Detected convergent synthesis with {complex_fragments} complex fragments (>8 atoms each) at depth {depth}"
                        )
                except Exception as e:
                    print(f"Error processing reaction: {e}")

        # Continue traversing
        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    # Start traversal
    dfs_traverse(route)
    return is_convergent
