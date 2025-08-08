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
    This function detects a synthetic strategy involving late-stage amide formation
    as the final step in a linear synthesis.
    """
    amide_formation_at_final_step = False
    is_linear_synthesis = True
    reaction_count = 0

    def dfs_traverse(node, depth=0, branch_id=0):
        nonlocal amide_formation_at_final_step, is_linear_synthesis, reaction_count

        print(f"Traversing node at depth {depth}, type: {node['type']}, branch: {branch_id}")

        # Check if this is a reaction node
        if node["type"] == "reaction" and "metadata" in node and "rsmi" in node["metadata"]:
            reaction_count += 1
            rsmi = node["metadata"]["rsmi"]
            reactants_smiles = rsmi.split(">")[0].split(".")
            product_smiles = rsmi.split(">")[-1]

            print(f"Checking reaction at depth {depth}: {rsmi}")

            # In retrosynthesis, depth 1 is the final step in forward synthesis
            is_final_step = depth == 1

            if is_final_step:
                print(f"This is the final step reaction (depth {depth})")

                # Check for amide formation reactions using the checker function
                amide_formation_reactions = [
                    "Acylation of Nitrogen Nucleophiles by Acyl/Thioacyl/Carbamoyl Halides and Analogs_N",
                    "Acylation of Nitrogen Nucleophiles by Carboxylic Acids",
                    "Acyl chloride with primary amine to amide (Schotten-Baumann)",
                    "Acyl chloride with secondary amine to amide",
                    "Carboxylic acid with primary amine to amide",
                    "Ester with primary amine to amide",
                    "Ester with secondary amine to amide",
                    "Ester with ammonia to amide",
                    "Acyl chloride with ammonia to amide",
                    "Schotten-Baumann_amide",
                ]

                is_amide_formation = False
                for rxn_type in amide_formation_reactions:
                    if checker.check_reaction(rxn_type, rsmi):
                        print(f"Detected amide formation reaction: {rxn_type}")
                        is_amide_formation = True
                        break

                if not is_amide_formation:
                    print(
                        "No specific amide formation reaction detected, checking functional groups..."
                    )
                    # Fallback to functional group checking
                    has_acid = any(checker.check_fg("Carboxylic acid", r) for r in reactants_smiles)
                    has_acyl_halide = any(
                        checker.check_fg("Acyl halide", r) for r in reactants_smiles
                    )
                    has_ester = any(checker.check_fg("Ester", r) for r in reactants_smiles)

                    has_primary_amine = any(
                        checker.check_fg("Primary amine", r) for r in reactants_smiles
                    )
                    has_secondary_amine = any(
                        checker.check_fg("Secondary amine", r) for r in reactants_smiles
                    )
                    has_aniline = any(checker.check_fg("Aniline", r) for r in reactants_smiles)

                    # Check if product has amide but reactants don't
                    has_primary_amide_product = checker.check_fg("Primary amide", product_smiles)
                    has_secondary_amide_product = checker.check_fg(
                        "Secondary amide", product_smiles
                    )
                    has_tertiary_amide_product = checker.check_fg("Tertiary amide", product_smiles)
                    has_amide_product = (
                        has_primary_amide_product
                        or has_secondary_amide_product
                        or has_tertiary_amide_product
                    )

                    has_amide_reactants = any(
                        checker.check_fg("Primary amide", r)
                        or checker.check_fg("Secondary amide", r)
                        or checker.check_fg("Tertiary amide", r)
                        for r in reactants_smiles
                    )

                    print(
                        f"FG analysis - Acid: {has_acid}, Acyl halide: {has_acyl_halide}, Ester: {has_ester}"
                    )
                    print(
                        f"FG analysis - Primary amine: {has_primary_amine}, Secondary amine: {has_secondary_amine}, Aniline: {has_aniline}"
                    )
                    print(
                        f"FG analysis - Product has amide: {has_amide_product}, Reactants have amide: {has_amide_reactants}"
                    )

                    # Check if we have the right reactants and product for amide formation
                    # and that the amide is actually formed in this step
                    if (
                        (has_acid or has_acyl_halide or has_ester)
                        and (has_primary_amine or has_secondary_amine or has_aniline)
                        and has_amide_product
                        and not has_amide_reactants
                    ):
                        print(
                            "Detected amide formation at final step through functional group analysis"
                        )
                        is_amide_formation = True

                if is_amide_formation:
                    print("Confirmed amide formation at final step")
                    amide_formation_at_final_step = True

        # Check for branching in the synthesis route
        if node["type"] == "mol" and len(node.get("children", [])) > 1:
            # More than one reaction from this molecule indicates branching
            print(f"Detected branching at depth {depth}")
            is_linear_synthesis = False

        # Continue traversing
        child_count = 0
        for child in node.get("children", []):
            child_count += 1
            # Pass a new branch_id for each child of a molecule with multiple children
            new_branch_id = branch_id
            if node["type"] == "mol" and child_count > 1:
                new_branch_id = branch_id + 1
            dfs_traverse(child, depth + 1, new_branch_id)

    # Start traversal from the root
    dfs_traverse(route)

    print(f"Final result: {amide_formation_at_final_step and is_linear_synthesis}")
    return amide_formation_at_final_step and is_linear_synthesis
