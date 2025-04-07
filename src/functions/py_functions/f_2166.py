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
    Detects if the synthesis route follows a semi-convergent approach with multiple
    fragment couplings rather than linear chain extension.
    """
    # Count reactions with multiple complex reactants
    multi_fragment_reactions = 0
    total_reactions = 0

    def dfs_traverse(node):
        nonlocal multi_fragment_reactions, total_reactions

        if (
            node["type"] == "reaction"
            and "metadata" in node
            and "rsmi" in node["metadata"]
        ):
            total_reactions += 1
            rsmi = node["metadata"]["rsmi"]
            reactants = rsmi.split(">")[0].split(".")

            # Count complex reactants (those with multiple rings or heterocycles)
            complex_reactants = 0
            for reactant in reactants:
                reactant_mol = Chem.MolFromSmiles(reactant)
                if not reactant_mol:
                    continue

                # Check if reactant has rings
                ring_count = Chem.GetSSSR(reactant_mol)
                if len(ring_count) >= 1:
                    complex_reactants += 1

            if complex_reactants >= 2:
                multi_fragment_reactions += 1
                print(
                    f"Found multi-fragment coupling at depth: {node.get('depth', 'unknown')}"
                )

        for child in node.get("children", []):
            dfs_traverse(child)

    dfs_traverse(route)

    # If at least 1/3 of reactions involve multiple complex fragments, consider it semi-convergent
    return total_reactions > 0 and (multi_fragment_reactions / total_reactions) >= 0.33
