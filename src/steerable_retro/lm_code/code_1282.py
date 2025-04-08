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

from steerable_retro.utils import check, fuzzy_dict
from steerable_retro.utils.check import Check

root_data = "/home/andres/Documents/steerable_retro/data"

fg_args = {
    "file_path": f"{root_data}/patterns/functional_groups.json",
    "value_field": "pattern",
    "key_field": "name",
}
reaction_class_args = {
    "file_path": f"{root_data}/patterns/smirks.json",
    "value_field": "smirks",
    "key_field": "name",
}
ring_smiles_args = {
    "file_path": f"{root_data}/patterns/chemical_rings_smiles.json",
    "value_field": "smiles",
    "key_field": "name",
}
functional_groups = fuzzy_dict.FuzzyDict.from_json(**fg_args)
reaction_classes = fuzzy_dict.FuzzyDict.from_json(**reaction_class_args)
ring_smiles = fuzzy_dict.FuzzyDict.from_json(**ring_smiles_args)

checker = check.Check(
    fg_dict=functional_groups, reaction_dict=reaction_classes, ring_dict=ring_smiles
)


def main(route):
    """
    Detects a functional group transformation sequence from trifluoroacetyl to carboxylic acid to hydrazide.
    In retrosynthetic traversal, we'll detect hydrazide first, then carboxylic acid, then trifluoroacetyl.
    """
    # Initialize tracking variables
    transformations = []

    def dfs_traverse(node, depth=0):
        nonlocal transformations

        if node["type"] == "reaction":
            if "metadata" in node and "rsmi" in node["metadata"]:
                rsmi = node["metadata"]["rsmi"]
                reactants = rsmi.split(">")[0].split(".")
                product = rsmi.split(">")[-1]

                # Check for hydrazide to carboxylic acid transformation (retrosynthetic)
                if checker.check_fg("Acylhydrazine", product) and any(
                    checker.check_fg("Carboxylic acid", r) for r in reactants
                ):
                    print(
                        f"Depth {depth}: Detected hydrazide to carboxylic acid transformation (retrosynthetic)"
                    )
                    transformations.append(("hydrazide", "carboxylic_acid", depth))

                # Check for carboxylic acid to trifluoroacetyl transformation (retrosynthetic)
                if checker.check_fg("Carboxylic acid", product) and any(
                    checker.check_fg("Trifluoro group", r) for r in reactants
                ):
                    print(
                        f"Depth {depth}: Detected carboxylic acid to trifluoroacetyl transformation (retrosynthetic)"
                    )
                    transformations.append(("carboxylic_acid", "trifluoroacetyl", depth))

        # Traverse children
        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    # Start traversal from the root
    dfs_traverse(route)

    # Sort transformations by depth to get the correct sequence
    transformations.sort(key=lambda x: x[2])

    # Check if we have both required transformations
    has_hydrazide_to_acid = False
    has_acid_to_trifluoro = False

    for from_fg, to_fg, _ in transformations:
        if from_fg == "hydrazide" and to_fg == "carboxylic_acid":
            has_hydrazide_to_acid = True
        if from_fg == "carboxylic_acid" and to_fg == "trifluoroacetyl":
            has_acid_to_trifluoro = True

    # Check if transformations are in the correct order
    correct_sequence = False
    if has_hydrazide_to_acid and has_acid_to_trifluoro:
        # Get indices to check sequence
        hydrazide_to_acid_idx = next(
            (
                i
                for i, (from_fg, to_fg, _) in enumerate(transformations)
                if from_fg == "hydrazide" and to_fg == "carboxylic_acid"
            ),
            -1,
        )
        acid_to_trifluoro_idx = next(
            (
                i
                for i, (from_fg, to_fg, _) in enumerate(transformations)
                if from_fg == "carboxylic_acid" and to_fg == "trifluoroacetyl"
            ),
            -1,
        )

        # In retrosynthetic order, hydrazide should be found before carboxylic acid
        if (
            hydrazide_to_acid_idx != -1
            and acid_to_trifluoro_idx != -1
            and hydrazide_to_acid_idx < acid_to_trifluoro_idx
        ):
            correct_sequence = True

    print(f"Transformations: {transformations}")
    print(f"Trifluoroacetyl → carboxylic acid → hydrazide sequence detected: {correct_sequence}")
    return correct_sequence
