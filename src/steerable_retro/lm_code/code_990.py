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
    This function detects if the final step (depth 0) is an amide formation.
    """
    late_stage_amide_formation_detected = False

    def dfs_traverse(node, depth=0):
        nonlocal late_stage_amide_formation_detected

        print(f"Traversing node type: {node['type']} at depth: {depth}")

        if node["type"] == "reaction" and depth <= 1:  # Check depth 0 and 1 for final reaction
            print(f"Examining reaction at depth {depth}")
            try:
                rsmi = node["metadata"]["rsmi"]
                reactants_smiles = rsmi.split(">")[0].split(".")
                product_smiles = rsmi.split(">")[-1]

                print(f"Reaction SMILES: {rsmi}")
                print(f"Reactants: {reactants_smiles}")
                print(f"Product SMILES: {product_smiles}")

                # Check if the reaction is an amide formation using the checker function
                amide_formation_reactions = [
                    "Acylation of Nitrogen Nucleophiles by Acyl/Thioacyl/Carbamoyl Halides and Analogs_N",
                    "Acylation of Nitrogen Nucleophiles by Acyl/Thioacyl/Carbamoyl Halides and Analogs_OS",
                    "Acylation of Nitrogen Nucleophiles by Carboxylic Acids",
                    "Acyl chloride with ammonia to amide",
                    "Acyl chloride with primary amine to amide (Schotten-Baumann)",
                    "Acyl chloride with secondary amine to amide",
                    "Carboxylic acid with primary amine to amide",
                    "Ester with ammonia to amide",
                    "Ester with primary amine to amide",
                    "Ester with secondary amine to amide",
                    "Schotten-Baumann_amide",
                ]

                is_amide_formation = any(
                    checker.check_reaction(rxn_type, rsmi) for rxn_type in amide_formation_reactions
                )

                if is_amide_formation:
                    print("Detected amide formation reaction via reaction type check")
                    late_stage_amide_formation_detected = True
                    return

                # If reaction type check fails, check for functional group patterns
                has_carboxylic_acid = any(
                    checker.check_fg("Carboxylic acid", reactant) for reactant in reactants_smiles
                )
                has_acyl_halide = any(
                    checker.check_fg("Acyl halide", reactant) for reactant in reactants_smiles
                )
                has_ester = any(
                    checker.check_fg("Ester", reactant) for reactant in reactants_smiles
                )

                has_primary_amine = any(
                    checker.check_fg("Primary amine", reactant) for reactant in reactants_smiles
                )
                has_secondary_amine = any(
                    checker.check_fg("Secondary amine", reactant) for reactant in reactants_smiles
                )
                has_ammonia = any(reactant.strip() == "N" for reactant in reactants_smiles)

                # Check if amide is in product but not in reactants (i.e., it was formed)
                has_amide_in_product = (
                    checker.check_fg("Primary amide", product_smiles)
                    or checker.check_fg("Secondary amide", product_smiles)
                    or checker.check_fg("Tertiary amide", product_smiles)
                )

                has_amide_in_reactants = any(
                    checker.check_fg("Primary amide", reactant)
                    or checker.check_fg("Secondary amide", reactant)
                    or checker.check_fg("Tertiary amide", reactant)
                    for reactant in reactants_smiles
                )

                has_new_amide = has_amide_in_product and not has_amide_in_reactants

                print(f"Has carboxylic acid: {has_carboxylic_acid}")
                print(f"Has acyl halide: {has_acyl_halide}")
                print(f"Has ester: {has_ester}")
                print(f"Has primary amine: {has_primary_amine}")
                print(f"Has secondary amine: {has_secondary_amine}")
                print(f"Has ammonia: {has_ammonia}")
                print(f"Has amide in product: {has_amide_in_product}")
                print(f"Has amide in reactants: {has_amide_in_reactants}")
                print(f"Has new amide formed: {has_new_amide}")

                # Check if we have the necessary components for amide formation
                if (
                    (has_carboxylic_acid or has_acyl_halide or has_ester)
                    and (has_primary_amine or has_secondary_amine or has_ammonia)
                    and has_new_amide
                ):
                    print("Detected late-stage amide formation based on functional groups")
                    late_stage_amide_formation_detected = True
            except Exception as e:
                print(f"Error analyzing reaction: {e}")

        # Traverse children with incremented depth
        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    dfs_traverse(route)
    return late_stage_amide_formation_detected
