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
    This function detects if the synthetic route follows a linear fragment assembly strategy
    rather than a convergent approach.
    """
    # Track number of reactions with multiple reactants (potential convergent steps)
    convergent_steps = 0
    total_steps = 0

    def is_reagent(smiles):
        """Helper function to identify common reagents that shouldn't count as complex fragments"""
        mol = Chem.MolFromSmiles(smiles)
        if not mol:
            return False

        # Small molecules are likely reagents
        if mol.GetNumHeavyAtoms() <= 6:
            return True

        # Check for common reagent functional groups
        common_reagents = [
            "Acyl halide",
            "Triflate",
            "Mesylate",
            "Tosylate",
            "Magnesium halide",
            "Zinc halide",
            "Primary halide",
            "Secondary halide",
            "Tertiary halide",
            "Boronic acid",
            "Boronic ester",
            "Sulfonyl halide",
        ]

        for reagent in common_reagents:
            if checker.check_fg(reagent, smiles):
                return True

        return False

    def is_convergent_reaction(rsmi):
        """Check if the reaction type is inherently convergent"""
        convergent_reaction_types = [
            "Suzuki coupling with boronic acids",
            "Suzuki coupling with boronic esters",
            "Buchwald-Hartwig/Ullmann-Goldberg/N-arylation primary amine",
            "Buchwald-Hartwig/Ullmann-Goldberg/N-arylation secondary amine",
            "Stille reaction_aryl",
            "Negishi coupling",
            "Heck terminal vinyl",
            "Sonogashira alkyne_aryl halide",
        ]

        for rxn_type in convergent_reaction_types:
            if checker.check_reaction(rxn_type, rsmi):
                return True

        return False

    def dfs_traverse(node):
        nonlocal convergent_steps, total_steps

        if node["type"] == "reaction":
            try:
                rsmi = node["metadata"]["rsmi"]
                reactants_smiles = rsmi.split(">")[0].split(".")

                total_steps += 1

                # Check if the reaction type is inherently convergent
                if is_convergent_reaction(rsmi):
                    convergent_steps += 1
                    print(f"Detected inherently convergent reaction type: {rsmi}")
                # If a reaction has more than one reactant, it might be a convergent step
                elif len(reactants_smiles) > 1:
                    # Check if both reactants are complex (not just reagents)
                    complex_reactants = 0
                    for reactant_smiles in reactants_smiles:
                        if not is_reagent(reactant_smiles):
                            complex_reactants += 1

                    if complex_reactants > 1:
                        convergent_steps += 1
                        print(f"Detected potential convergent step: {rsmi}")
            except Exception as e:
                print(f"Error processing reaction node: {e}")

        for child in node.get("children", []):
            dfs_traverse(child)

    dfs_traverse(route)

    # Handle edge case of no reaction steps
    if total_steps == 0:
        print("No reaction steps found")
        return False

    # If 25% or less of steps are convergent, consider it a linear strategy
    is_linear = (convergent_steps / total_steps) <= 0.25

    print(f"Total steps: {total_steps}, Convergent steps: {convergent_steps}")
    print(f"Is linear fragment assembly: {is_linear}")

    return is_linear
