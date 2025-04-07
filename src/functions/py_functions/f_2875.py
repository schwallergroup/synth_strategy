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

fg_args = {
    "file_path": "/home/dparm/steerable_retro/data/patterns/functional_groups.json",
    "value_field": "pattern",
    "key_field": "name",
}
reaction_class_args = {
    "file_path": "/home/dparm/steerable_retro/data/patterns/smirks.json",
    "value_field": "smirks",
    "key_field": "name",
}
ring_smiles_args = {
    "file_path": "/home/dparm/steerable_retro/data/patterns/chemical_rings_smiles.json",
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
    Detects if the synthesis includes a Grignard addition to form a C-C bond.
    """
    grignard_detected = False

    def dfs_traverse(node):
        nonlocal grignard_detected

        if node["type"] == "reaction" and not grignard_detected:
            if "rsmi" in node.get("metadata", {}):
                rsmi = node["metadata"]["rsmi"]

                # Check for Grignard addition reactions using the checker function
                if (
                    checker.check_reaction("Grignard from aldehyde to alcohol", rsmi)
                    or checker.check_reaction("Grignard from ketone to alcohol", rsmi)
                    or checker.check_reaction(
                        "Grignard with CO2 to carboxylic acid", rsmi
                    )
                    or checker.check_reaction("Grignard from nitrile to ketone", rsmi)
                    or checker.check_reaction(
                        "Olefination of ketones with Grignard reagents", rsmi
                    )
                    or checker.check_reaction(
                        "Olefination of aldehydes with Grignard reagents", rsmi
                    )
                    or checker.check_reaction("Formation of Grignard reagents", rsmi)
                ):
                    print(f"Grignard addition detected in reaction: {rsmi}")
                    grignard_detected = True

                # Additional check for reactions involving Grignard reagents
                if not grignard_detected:
                    try:
                        reactants = rsmi.split(">")[0].split(".")
                        product = rsmi.split(">")[-1]

                        # Check if any reactant contains a magnesium halide (Grignard reagent)
                        for reactant in reactants:
                            if checker.check_fg("Magnesium halide", reactant):
                                print(
                                    f"Grignard reagent detected in reactant: {reactant}"
                                )
                                grignard_detected = True
                                break
                    except Exception as e:
                        print(f"Error analyzing reaction components: {e}")

        # Continue traversal even if we've found a Grignard reaction
        for child in node.get("children", []):
            dfs_traverse(child)

    dfs_traverse(route)
    print(f"Grignard addition strategy: {grignard_detected}")
    return grignard_detected
