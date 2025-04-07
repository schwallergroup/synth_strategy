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
    This function detects if the synthesis follows a linear strategy (no convergent steps with multiple complex fragments until final step).

    In retrosynthetic analysis:
    - Depth 0: Target molecule
    - Depth 1: Final synthetic step (convergent steps allowed here)
    - Depth >1: Earlier steps (should be linear)
    """
    is_linear = True

    # Handle empty routes
    if not route or "children" not in route or not route.get("children"):
        print("Route has no children, considering it linear by default")
        return True

    def dfs_traverse(node, depth=0):
        nonlocal is_linear

        print(f"Traversing node at depth {depth}, type: {node['type']}")

        if node["type"] == "reaction":
            if "rsmi" in node.get("metadata", {}):
                rsmi = node["metadata"]["rsmi"]
                reactants_smiles = rsmi.split(">")[0].split(".")

                # Check if this is a convergent step (multiple complex reactants)
                complex_reactants = 0
                for reactant in reactants_smiles:
                    reactant_mol = Chem.MolFromSmiles(reactant)
                    if reactant_mol:
                        # Define "complex" as having more than 7 heavy atoms
                        if reactant_mol.GetNumHeavyAtoms() > 7:
                            complex_reactants += 1
                            print(
                                f"Complex reactant found at depth {depth}: {reactant} with {reactant_mol.GetNumHeavyAtoms()} heavy atoms"
                            )

                # If more than one complex reactant, it's a convergent step
                if complex_reactants > 1:
                    print(
                        f"Convergent step detected at depth {depth} with {complex_reactants} complex reactants"
                    )
                    # If not the final step (depth > 1), mark as non-linear
                    # In retrosynthetic analysis, depth 1 is the final synthetic step
                    if depth > 1:
                        is_linear = False
                        print(
                            f"Non-linear synthesis: convergent step at depth {depth} (not the final step)"
                        )

        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    dfs_traverse(route)
    print(
        f"Final determination: Synthesis is {'linear' if is_linear else 'non-linear'}"
    )
    return is_linear
