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
    This function detects if borylation is used to prepare a coupling partner
    before a Suzuki coupling reaction.
    """
    borylation_detected = False
    suzuki_depth = None
    borylation_depth = None

    def dfs_traverse(node, depth=0):
        nonlocal borylation_detected, suzuki_depth, borylation_depth

        if (
            node["type"] == "reaction"
            and "metadata" in node
            and "rsmi" in node["metadata"]
        ):
            rsmi = node["metadata"]["rsmi"]
            reactants = rsmi.split(">")[0].split(".")
            product = rsmi.split(">")[-1]

            # Check for Suzuki coupling
            if len(reactants) > 1:
                boronic_pattern = Chem.MolFromSmarts("[#6]-[#5;X3,X4]")
                halide_pattern = Chem.MolFromSmarts("[#6]:[#6]-[#35,#53,#17]")

                has_boronic = False
                has_halide = False

                for reactant in reactants:
                    try:
                        mol = Chem.MolFromSmiles(reactant)
                        if mol:
                            if mol.HasSubstructMatch(boronic_pattern):
                                has_boronic = True
                            if mol.HasSubstructMatch(halide_pattern):
                                has_halide = True
                    except:
                        continue

                if has_boronic and has_halide:
                    suzuki_depth = depth

            # Check for borylation reaction
            if len(reactants) >= 1:
                # Look for conversion of aryl halide to boronic ester
                try:
                    for reactant in reactants:
                        mol_reactant = Chem.MolFromSmiles(reactant)
                        if mol_reactant and mol_reactant.HasSubstructMatch(
                            Chem.MolFromSmarts("[#6]:[#6]-[#35,#53,#17]")
                        ):
                            mol_product = Chem.MolFromSmiles(product)
                            if mol_product and mol_product.HasSubstructMatch(
                                Chem.MolFromSmarts("[#6]:[#6]-[#5]")
                            ):
                                print(
                                    "Borylation detected in reaction with RSMI:", rsmi
                                )
                                borylation_depth = depth
                except:
                    pass

        # Continue traversing
        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    # Start traversal
    dfs_traverse(route)

    # Check if borylation occurs before Suzuki coupling
    if suzuki_depth is not None and borylation_depth is not None:
        if (
            borylation_depth > suzuki_depth
        ):  # Higher depth means earlier in synthesis (retrosynthetic direction)
            borylation_detected = True
            print(
                f"Borylation (depth {borylation_depth}) occurs before Suzuki coupling (depth {suzuki_depth})"
            )

    return borylation_detected
