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
    This function detects a linear synthesis with late-stage C-C bond formation via alkylation,
    preceded by functional group interconversions (ester to alcohol to halide).
    """
    # Track key transformations and their sequence
    transformations = []

    def dfs_traverse(node, depth=0):
        nonlocal transformations

        if node["type"] == "reaction":
            # Extract reaction SMILES
            rsmi = node["metadata"]["rsmi"]
            reactants_smiles = rsmi.split(">")[0].split(".")
            product_smiles = rsmi.split(">")[-1]

            print(f"Analyzing reaction at depth {depth}: {rsmi}")

            # Check for ester reduction to alcohol
            ester_reduction_reactions = [
                "Reduction of ester to primary alcohol",
                "Reduction of carboxylic acid to primary alcohol",
            ]

            if any(checker.check_reaction(rxn, rsmi) for rxn in ester_reduction_reactions):
                print(f"Detected ester reduction at depth {depth}")
                transformations.append(("ester_reduction", depth))
            # Additional check for ester to alcohol conversion
            elif any(checker.check_fg("Ester", r) for r in reactants_smiles) and checker.check_fg(
                "Primary alcohol", product_smiles
            ):
                print(f"Detected potential ester reduction at depth {depth}")
                transformations.append(("ester_reduction", depth))

            # Check for alcohol to halide conversion
            halide_reactions = [
                "Alcohol to chloride_SOCl2",
                "Alcohol to chloride_PCl5_ortho",
                "Alcohol to chloride_POCl3",
                "Alcohol to chloride_HCl",
                "Alcohol to bromide",
                "Alcohol to iodide",
                "Appel reaction",
            ]

            if any(checker.check_reaction(rxn, rsmi) for rxn in halide_reactions):
                print(f"Detected alcohol to halide conversion at depth {depth}")
                transformations.append(("alcohol_to_halide", depth))
            # Additional check for alcohol to halide conversion
            elif any(checker.check_fg("Primary alcohol", r) for r in reactants_smiles) and any(
                checker.check_fg(halide, product_smiles)
                for halide in ["Primary halide", "Secondary halide", "Tertiary halide"]
            ):
                print(f"Detected potential alcohol to halide conversion at depth {depth}")
                transformations.append(("alcohol_to_halide", depth))

            # Check for alkylation (C-C bond formation)
            alkylation_reactions = [
                "Alkylation of amines",
                "N-alkylation of primary amines with alkyl halides",
                "N-alkylation of secondary amines with alkyl halides",
                "S-alkylation of thiols",
                "S-alkylation of thiols with alcohols",
                "Friedel-Crafts alkylation",
                "Friedel-Crafts alkylation with halide",
            ]

            if any(checker.check_reaction(rxn, rsmi) for rxn in alkylation_reactions):
                print(f"Detected alkylation reaction at depth {depth}")
                transformations.append(("alkylation", depth))
            # Check for C-C bond formation with halide
            elif any(
                any(
                    checker.check_fg(halide, r)
                    for halide in ["Primary halide", "Secondary halide", "Tertiary halide"]
                )
                for r in reactants_smiles
            ):
                # This is a simplified check for C-C bond formation
                # In a real scenario, we would need more sophisticated analysis
                print(f"Detected potential C-C bond formation with halide at depth {depth}")
                transformations.append(("alkylation", depth))

        # Traverse children
        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    # Start traversal
    dfs_traverse(route)

    # Sort transformations by depth to get the sequence
    transformations.sort(key=lambda x: x[1])
    reaction_sequence = [t[0] for t in transformations]

    print(f"Detected transformations: {transformations}")
    print(f"Reaction sequence: {reaction_sequence}")

    # Check if we have all three required transformations
    has_ester_reduction = "ester_reduction" in reaction_sequence
    has_alcohol_to_halide = "alcohol_to_halide" in reaction_sequence
    has_alkylation = "alkylation" in reaction_sequence

    # Check if they appear in the correct order
    correct_order = False
    if has_ester_reduction and has_alcohol_to_halide and has_alkylation:
        ester_idx = reaction_sequence.index("ester_reduction")
        halide_idx = reaction_sequence.index("alcohol_to_halide")
        alkylation_idx = reaction_sequence.index("alkylation")

        # In retrosynthetic traversal, higher depth means earlier in synthesis
        # So alkylation (late stage) should have lowest depth, followed by halide formation, then ester reduction
        correct_order = alkylation_idx < halide_idx < ester_idx

    all_present = has_ester_reduction and has_alcohol_to_halide and has_alkylation

    print(f"All transformations present: {all_present}")
    print(f"Correct order: {correct_order}")

    return all_present and correct_order
