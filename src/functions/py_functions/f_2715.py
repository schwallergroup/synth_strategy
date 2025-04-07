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
    Detects if the synthesis involves a late-stage N-oxidation of a pyridine or similar heterocycle.
    Late stage is defined as occurring in the first half of the synthesis (low depth).
    """
    n_oxidation_detected = False
    max_depth = 0
    n_oxidation_depth = None

    def dfs_traverse(node, depth=0):
        nonlocal n_oxidation_detected, max_depth, n_oxidation_depth

        max_depth = max(max_depth, depth)

        if node["type"] == "reaction":
            # Check if this is an N-oxidation reaction
            if "rsmi" in node.get("metadata", {}):
                rsmi = node["metadata"]["rsmi"]
                reactants = rsmi.split(">")[0]
                product = rsmi.split(">")[-1]

                # Check if product has N-oxide but reactant doesn't
                product_mol = Chem.MolFromSmiles(product)
                reactant_mols = [Chem.MolFromSmiles(r) for r in reactants.split(".")]

                if product_mol and all(reactant_mol for reactant_mol in reactant_mols):
                    n_oxide_pattern = Chem.MolFromSmarts("[n+][O-]")
                    if product_mol.HasSubstructMatch(n_oxide_pattern):
                        # Check if any reactant has N-oxide
                        reactant_has_n_oxide = any(
                            r.HasSubstructMatch(n_oxide_pattern)
                            for r in reactant_mols
                            if r
                        )

                        if not reactant_has_n_oxide:
                            n_oxidation_detected = True
                            n_oxidation_depth = depth
                            print(f"N-oxidation detected at depth {depth}")

        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    dfs_traverse(route)

    # Consider it late-stage if it occurs in the first half of the synthesis
    is_late_stage = n_oxidation_depth is not None and n_oxidation_depth <= max_depth / 2

    if is_late_stage:
        print(
            f"Late-stage N-oxidation detected at depth {n_oxidation_depth} (max depth: {max_depth})"
        )

    return n_oxidation_detected and is_late_stage
