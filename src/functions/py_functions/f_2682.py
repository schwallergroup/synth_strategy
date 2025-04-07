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
    This function detects convergent synthesis where two complex fragments are joined in a late stage.
    """
    convergent_synthesis = False

    def dfs_traverse(node):
        nonlocal convergent_synthesis

        # Check only reactions at depth 0 or 1 (late stage)
        if (
            node["type"] == "reaction"
            and "metadata" in node
            and "rsmi" in node["metadata"]
        ):
            depth = node["metadata"].get("depth", -1)
            if depth <= 1:  # Late stage reaction
                rsmi = node["metadata"]["rsmi"]
                reactants = rsmi.split(">")[0].split(".")

                # Check if we have at least 2 complex reactants
                if len(reactants) >= 2:
                    complex_reactants = 0

                    for reactant in reactants:
                        reactant_mol = Chem.MolFromSmiles(reactant)
                        if reactant_mol:
                            # Define complexity as having at least 10 atoms and 1 ring
                            atom_count = reactant_mol.GetNumAtoms()
                            ring_count = reactant_mol.GetRingInfo().NumRings()

                            if atom_count >= 10 and ring_count >= 1:
                                complex_reactants += 1

                    if complex_reactants >= 2:
                        convergent_synthesis = True
                        print("Convergent synthesis with complex fragments detected")

        # Traverse children
        for child in node.get("children", []):
            dfs_traverse(child)

    # Start traversal from the root
    dfs_traverse(route)

    return convergent_synthesis
