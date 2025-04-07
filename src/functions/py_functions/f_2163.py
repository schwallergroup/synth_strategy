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
    Detects if the synthesis route involves formation of a sulfonamide bond in the late stage
    (depth 0 or 1).
    """
    found_late_stage_sulfonamide = False

    def dfs_traverse(node):
        nonlocal found_late_stage_sulfonamide

        if (
            node["type"] == "reaction"
            and "metadata" in node
            and "rsmi" in node["metadata"]
        ):
            depth = node.get("depth", -1)
            if depth <= 1:  # Late stage (depth 0 or 1)
                rsmi = node["metadata"]["rsmi"]
                reactants = rsmi.split(">")[0].split(".")
                product = rsmi.split(">")[-1]

                # Check if product has sulfonamide but reactants don't
                sulfonamide_pattern = Chem.MolFromSmarts("[#16](=[O])(=[O])[#7]")

                product_mol = Chem.MolFromSmiles(product)
                if not product_mol:
                    return

                if product_mol.HasSubstructMatch(sulfonamide_pattern):
                    # Check if any reactant has this pattern
                    has_pattern_in_reactants = False
                    for reactant in reactants:
                        reactant_mol = Chem.MolFromSmiles(reactant)
                        if not reactant_mol:
                            continue
                        if reactant_mol.HasSubstructMatch(sulfonamide_pattern):
                            has_pattern_in_reactants = True
                            break

                    if not has_pattern_in_reactants:
                        found_late_stage_sulfonamide = True
                        print(
                            f"Found late-stage sulfonamide formation at depth: {depth}"
                        )

        for child in node.get("children", []):
            dfs_traverse(child)

    dfs_traverse(route)

    return found_late_stage_sulfonamide
