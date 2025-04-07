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
    Detects a linear synthesis strategy (as opposed to convergent) where most steps
    involve only one or two reactants and build complexity sequentially.
    """
    # Track reactant counts at each step
    reactant_counts = []
    total_steps = 0

    def dfs_traverse(node):
        nonlocal total_steps

        if node["type"] == "reaction":
            # Extract reactants
            rsmi = node["metadata"]["rsmi"]
            reactants = rsmi.split(">")[0].split(".")

            # Count non-trivial reactants (more than 3 atoms)
            significant_reactants = 0
            for r in reactants:
                try:
                    mol = Chem.MolFromSmiles(r)
                    if mol and mol.GetNumAtoms() > 3:
                        significant_reactants += 1
                except:
                    pass

            reactant_counts.append(significant_reactants)
            total_steps += 1

        # Process children
        for child in node.get("children", []):
            dfs_traverse(child)

    # Start traversal
    dfs_traverse(route)

    # Analyze reactant counts
    if not reactant_counts:
        return False

    # Calculate metrics
    avg_reactants = sum(reactant_counts) / len(reactant_counts)
    max_reactants = max(reactant_counts)
    steps_with_one_or_two = sum(1 for count in reactant_counts if count <= 2)

    # Linear synthesis typically has most steps with 1-2 significant reactants
    is_linear = (
        avg_reactants < 2.5 and max_reactants <= 3 and steps_with_one_or_two >= 0.7 * total_steps
    )

    print(f"Linear synthesis strategy detection results:")
    print(f"  Reactant counts per step: {reactant_counts}")
    print(f"  Average reactants per step: {avg_reactants:.2f}")
    print(f"  Maximum reactants in any step: {max_reactants}")
    print(f"  Steps with 1-2 reactants: {steps_with_one_or_two}/{total_steps}")
    print(f"  Is linear synthesis: {is_linear}")

    return is_linear
