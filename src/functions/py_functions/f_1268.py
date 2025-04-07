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
    This function detects a strategy where a phenol group is preserved throughout most of the synthesis.
    """
    steps_with_phenol = 0
    total_steps = 0

    def dfs_traverse(node):
        nonlocal steps_with_phenol, total_steps

        if node["type"] == "reaction":
            total_steps += 1
            if "metadata" in node and "rsmi" in node["metadata"]:
                rsmi = node["metadata"]["rsmi"]
                reactants = rsmi.split(">")[0]

                reactant_mol = Chem.MolFromSmiles(reactants)

                if reactant_mol:
                    # Phenol pattern
                    phenol_pattern = Chem.MolFromSmarts("[OH]-[c]")

                    if reactant_mol.HasSubstructMatch(phenol_pattern):
                        steps_with_phenol += 1
                        print(f"Found phenol in step: {rsmi}")

        # Traverse children
        for child in node.get("children", []):
            dfs_traverse(child)

    # Start traversal
    dfs_traverse(route)

    # Check if phenol is preserved in most steps (>= 75% of steps)
    phenol_preservation_ratio = (
        steps_with_phenol / total_steps if total_steps > 0 else 0
    )
    is_phenol_preserved = phenol_preservation_ratio >= 0.75

    print(
        f"Steps with phenol: {steps_with_phenol}/{total_steps} ({phenol_preservation_ratio:.2f})"
    )
    print(f"Phenol preservation strategy: {is_phenol_preserved}")

    return is_phenol_preserved
