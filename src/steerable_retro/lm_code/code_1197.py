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
    Detects if the synthesis follows a linear strategy without convergent steps.
    This is indicated by each reaction having only one non-reagent reactant.

    A linear synthesis is characterized by:
    1. Each reaction node has exactly one non-starting material molecule child
    2. No convergent reaction types (e.g., coupling reactions)
    """
    is_linear = True

    def is_likely_reagent(smiles):
        """Helper function to identify common reagents/solvents/catalysts"""
        mol = Chem.MolFromSmiles(smiles)
        if not mol:
            return False

        # Small molecules are likely reagents
        if mol.GetNumHeavyAtoms() < 6:
            return True

        # Check for common solvents and reagents
        common_reagents = ["C1CCOC1", "O=C(O)O", "[K+]", "CC(=O)O", "CCO", "CCOC(=O)O"]
        if smiles in common_reagents:
            return True

        return False

    def is_heterocycle_formation(rsmi):
        """Check if the reaction is a heterocycle formation which appears convergent but is considered linear"""
        # Check for imidazole formation
        if checker.check_reaction("imidazole", rsmi) or checker.check_reaction("{imidazole}", rsmi):
            return True

        # Check for benzimidazole formation
        if (
            checker.check_reaction("benzimidazole formation from aldehyde", rsmi)
            or checker.check_reaction("benzimidazole formation from acyl halide", rsmi)
            or checker.check_reaction("benzimidazole formation from ester/carboxylic acid", rsmi)
            or checker.check_reaction("{benzimidazole_derivatives_carboxylic-acid/ester}", rsmi)
            or checker.check_reaction("{benzimidazole_derivatives_aldehyde}", rsmi)
        ):
            return True

        # Check for other heterocycle formations that are considered linear
        linear_heterocycles = [
            "Paal-Knorr pyrrole synthesis",
            "{Paal-Knorr pyrrole}",
            "{benzoxazole_arom-aldehyde}",
            "{benzoxazole_carboxylic-acid}",
            "{thiazole}",
            "{pyrazole}",
            "pyrazole formation",
        ]

        for rxn_type in linear_heterocycles:
            if checker.check_reaction(rxn_type, rsmi):
                return True

        return False

    def dfs_traverse(node, depth=0):
        nonlocal is_linear

        # If we've already determined it's not linear, no need to continue
        if not is_linear:
            return

        if node["type"] == "reaction":
            # Check reaction type for inherently convergent reactions
            if "metadata" in node and "rsmi" in node["metadata"]:
                rsmi = node["metadata"]["rsmi"]

                # Skip heterocycle formation reactions which appear convergent but are considered linear
                if is_heterocycle_formation(rsmi):
                    print(f"Detected heterocycle formation (considered linear): {rsmi}")
                    return

                # Check for coupling reactions which are inherently convergent
                convergent_reactions = [
                    "Suzuki coupling with boronic acids",
                    "Negishi coupling",
                    "Stille reaction_aryl",
                    "Heck terminal vinyl",
                    "Sonogashira alkyne_aryl halide",
                    "Buchwald-Hartwig",
                    "{Buchwald-Hartwig}",
                    "N-arylation (Buchwald-Hartwig/Ullmann-Goldberg)",
                    "Ugi reaction",
                ]

                for rxn_type in convergent_reactions:
                    if checker.check_reaction(rxn_type, rsmi):
                        print(f"Detected convergent reaction type: {rxn_type} in {rsmi}")
                        is_linear = False
                        return

                # Count non-starting material molecule children
                non_starting_material_count = 0
                for child in node.get("children", []):
                    if child["type"] == "mol" and not child.get("in_stock", False):
                        non_starting_material_count += 1

                # In a linear synthesis, each reaction should have exactly one
                # non-starting material molecule child
                if non_starting_material_count > 1:
                    print(
                        f"Detected convergent step with {non_starting_material_count} non-starting material reactants"
                    )
                    is_linear = False
                    return

                # Additional check: if there are multiple complex reactants in the RSMI
                reactants = rsmi.split(">")[0].split(".")
                complex_reactant_count = 0

                for reactant in reactants:
                    if is_likely_reagent(reactant):
                        continue

                    reactant_mol = Chem.MolFromSmiles(reactant)
                    if reactant_mol and reactant_mol.GetNumHeavyAtoms() > 12:  # Increased threshold
                        complex_reactant_count += 1

                # Check if this is a convergent step but exclude certain reaction types
                # that appear convergent but are considered linear in synthesis planning
                if complex_reactant_count > 1:
                    # Check if this is a common linear reaction type that might appear convergent
                    linear_reaction_types = [
                        "Acylation of Nitrogen Nucleophiles by Acyl/Thioacyl/Carbamoyl Halides and Analogs_N",
                        "Esterification of Carboxylic Acids",
                        "Acylation of primary amines",
                        "Acylation of secondary amines",
                        "Schotten-Baumann to ester",
                        "{Schotten-Baumann_amide}",
                    ]

                    if not any(checker.check_reaction(rxn, rsmi) for rxn in linear_reaction_types):
                        print(
                            f"Detected convergent step with multiple complex reactants in RSMI: {rsmi}"
                        )
                        is_linear = False
                        return

        # Traverse children
        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    # Start traversal
    dfs_traverse(route)
    return is_linear
