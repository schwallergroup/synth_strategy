#!/bin/python

"""LM-defined function for strategy description."""

from rdkit.Chem import AllChem, rdFMCS
import copy
from collections import deque
import rdkit.Chem as Chem
from rdkit.Chem import rdMolDescriptors
from rdkit.Chem import rdChemReactions
from rdkit.Chem import AllChem
from rdkit.Chem import rdFMCS
import rdkit.Chem.rdFMCS
from rdkit.Chem import AllChem, Descriptors, rdMolDescriptors
from rdkit import Chem
from rdkit.Chem import Descriptors
from rdkit.Chem import AllChem, rdMolDescriptors
from rdkit.Chem import AllChem, Descriptors, Lipinski
from rdkit.Chem import rdmolops
import re
from rdkit.Chem.Scaffolds import MurckoScaffold
from rdkit.Chem import AllChem, Descriptors
import traceback
import rdkit
from collections import Counter
from steerable_retro.utils.check import Check
from steerable_retro.utils import fuzzy_dict, check

root_data = "/home/dparm/steerable_retro/data"

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
    Detects a convergent synthesis strategy involving fragment coupling via amine formation.
    """
    # Track if we've found the required elements
    amine_formation_found = False

    def dfs_traverse(node):
        nonlocal amine_formation_found

        if node["type"] == "reaction":
            # Extract reactants and products
            try:
                rsmi = node["metadata"]["rsmi"]
                reactants_smiles = rsmi.split(">")[0].split(".")
                product_smiles = rsmi.split(">")[-1]

                # Check if this is a retrosynthetic disconnection of a secondary amine
                # (which means in forward synthesis, it's an amine formation)
                if len(reactants_smiles) >= 2:  # At least two fragments
                    # Check if product contains a secondary amine
                    if checker.check_fg("Secondary amine", product_smiles):
                        print(f"Found product with secondary amine: {product_smiles}")

                        # Check for various reactant patterns that could form secondary amines
                        has_aldehyde = any(
                            checker.check_fg("Aldehyde", smi) for smi in reactants_smiles
                        )
                        has_ketone = any(
                            checker.check_fg("Ketone", smi) for smi in reactants_smiles
                        )
                        has_primary_amine = any(
                            checker.check_fg("Primary amine", smi) for smi in reactants_smiles
                        )
                        has_primary_halide = any(
                            checker.check_fg("Primary halide", smi) for smi in reactants_smiles
                        )
                        has_secondary_halide = any(
                            checker.check_fg("Secondary halide", smi) for smi in reactants_smiles
                        )

                        print(
                            f"Has aldehyde: {has_aldehyde}, Has ketone: {has_ketone}, Has primary amine: {has_primary_amine}"
                        )
                        print(
                            f"Has primary halide: {has_primary_halide}, Has secondary halide: {has_secondary_halide}"
                        )

                        # Check for various amine formation reactions
                        is_reductive_amination_aldehyde = checker.check_reaction(
                            "Reductive amination with aldehyde", rsmi
                        )
                        is_reductive_amination_ketone = checker.check_reaction(
                            "Reductive amination with ketone", rsmi
                        )
                        is_addition_aldehyde = checker.check_reaction(
                            "Addition of primary amines to aldehydes/thiocarbonyls", rsmi
                        )
                        is_addition_ketone = checker.check_reaction(
                            "Addition of primary amines to ketones/thiocarbonyls", rsmi
                        )
                        is_n_alkylation = checker.check_reaction(
                            "N-alkylation of primary amines with alkyl halides", rsmi
                        )
                        is_buchwald_hartwig = checker.check_reaction(
                            "Buchwald-Hartwig/Ullmann-Goldberg/N-arylation primary amine", rsmi
                        )

                        print(f"Reductive amination (aldehyde): {is_reductive_amination_aldehyde}")
                        print(f"Reductive amination (ketone): {is_reductive_amination_ketone}")
                        print(f"Addition (aldehyde): {is_addition_aldehyde}")
                        print(f"Addition (ketone): {is_addition_ketone}")
                        print(f"N-alkylation: {is_n_alkylation}")
                        print(f"Buchwald-Hartwig: {is_buchwald_hartwig}")

                        # Check if any amine formation reaction is detected
                        if (
                            is_reductive_amination_aldehyde
                            or is_reductive_amination_ketone
                            or is_addition_aldehyde
                            or is_addition_ketone
                            or is_n_alkylation
                            or is_buchwald_hartwig
                        ):
                            amine_formation_found = True
                            print(f"Found amine formation via fragment coupling: {rsmi}")
                        # Fallback check based on reactant patterns if specific reaction not detected
                        elif ((has_aldehyde or has_ketone) and has_primary_amine) or (
                            has_primary_amine and (has_primary_halide or has_secondary_halide)
                        ):
                            # This is likely an amine formation reaction even if not specifically recognized
                            amine_formation_found = True
                            print(
                                f"Found likely amine formation based on reactant patterns: {rsmi}"
                            )
            except Exception as e:
                print(f"Error processing reaction node: {e}")

        # Traverse children
        for child in node.get("children", []):
            dfs_traverse(child)

    # Start traversal
    dfs_traverse(route)

    print(f"Convergent synthesis with amine formation strategy detected: {amine_formation_found}")
    return amine_formation_found
