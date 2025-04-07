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
    This function detects phenol alkylation strategy.
    Looks for transformation of phenol to phenolic ether.
    """
    phenol_alkylation_detected = False

    def dfs_traverse(node):
        nonlocal phenol_alkylation_detected

        if node["type"] == "reaction":
            if "rsmi" in node.get("metadata", {}):
                rsmi = node["metadata"]["rsmi"]
                reactants = rsmi.split(">")[0].split(".")
                product = rsmi.split(">")[-1]

                # Check for various phenol alkylation reactions
                if (
                    checker.check_reaction("Williamson Ether Synthesis", rsmi)
                    or checker.check_reaction("Mitsunobu aryl ether", rsmi)
                    or checker.check_reaction("Mitsunobu_phenole", rsmi)
                    or checker.check_reaction("{Williamson ether}", rsmi)
                ):
                    print(f"Detected potential phenol alkylation reaction: {rsmi}")

                    # Check for phenol in reactants
                    has_phenol = False
                    has_halide = False

                    for reactant in reactants:
                        if checker.check_fg("Phenol", reactant) or checker.check_fg(
                            "Aromatic alcohol", reactant
                        ):
                            print(f"Found phenol in reactant: {reactant}")
                            has_phenol = True

                        # Check for primary, secondary, or tertiary halide
                        if (
                            checker.check_fg("Primary halide", reactant)
                            or checker.check_fg("Secondary halide", reactant)
                            or checker.check_fg("Tertiary halide", reactant)
                        ):
                            print(f"Found halide in reactant: {reactant}")
                            has_halide = True

                    # Check for ether in product and no phenol in product
                    if (
                        has_phenol
                        and has_halide
                        and checker.check_fg("Ether", product)
                        and not checker.check_fg("Phenol", product)
                    ):
                        print(
                            f"Found ether in product and phenol was consumed: {product}"
                        )
                        phenol_alkylation_detected = True

                # Alternative check for phenol alkylation if not detected by reaction type
                if not phenol_alkylation_detected:
                    has_phenol = False
                    has_electrophile = False

                    for reactant in reactants:
                        if checker.check_fg("Phenol", reactant) or checker.check_fg(
                            "Aromatic alcohol", reactant
                        ):
                            print(f"Found phenol in reactant: {reactant}")
                            has_phenol = True

                        # Check for various electrophiles
                        if (
                            checker.check_fg("Primary halide", reactant)
                            or checker.check_fg("Secondary halide", reactant)
                            or checker.check_fg("Tertiary halide", reactant)
                            or checker.check_fg("Triflate", reactant)
                            or checker.check_fg("Mesylate", reactant)
                            or checker.check_fg("Tosylate", reactant)
                        ):
                            print(f"Found electrophile in reactant: {reactant}")
                            has_electrophile = True

                    # Check for ether formation and phenol consumption
                    if has_phenol and has_electrophile:
                        if checker.check_fg("Ether", product) and not checker.check_fg(
                            "Phenol", product
                        ):
                            print(
                                f"Detected phenol alkylation: phenol + electrophile → ether"
                            )
                            phenol_alkylation_detected = True

                # Check for general O-alkylation pattern
                if not phenol_alkylation_detected:
                    for reactant in reactants:
                        # If we have a phenol in reactants but not in product, and an ether in product
                        if (
                            checker.check_fg("Phenol", reactant)
                            and not checker.check_fg("Phenol", product)
                            and checker.check_fg("Ether", product)
                        ):
                            print(
                                f"Detected phenol transformation to ether: {reactant} → {product}"
                            )
                            phenol_alkylation_detected = True
                            break

        for child in node.get("children", []):
            dfs_traverse(child)

    dfs_traverse(route)
    return phenol_alkylation_detected
