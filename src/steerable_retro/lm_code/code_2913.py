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
    Detects a synthetic strategy involving late-stage amide coupling
    as the final step in the synthesis.
    """
    has_late_stage_amide = False

    def dfs_traverse(node, depth=0):
        nonlocal has_late_stage_amide

        print(f"Traversing node at depth {depth}, type: {node['type']}")

        if node["type"] == "reaction":
            # Check if this is the final step (depth 0 or 1)
            if depth <= 1:
                print(f"Checking potential late-stage reaction at depth {depth}")

                # Get reaction SMILES
                rsmi = node.get("metadata", {}).get("rsmi", "")
                if not rsmi:
                    print("No reaction SMILES found")
                    return

                print(f"Reaction SMILES: {rsmi}")

                # Extract reactants and product
                try:
                    reactants_part = rsmi.split(">")[0]
                    reactants_smiles = reactants_part.split(".")
                    product_smiles = rsmi.split(">")[-1]
                    print(f"Product: {product_smiles}")
                    print(f"Reactants: {', '.join(reactants_smiles)}")
                except Exception as e:
                    print(f"Error extracting reactants and product: {e}")
                    return

                # Check if product contains an amide
                has_amide_in_product = (
                    checker.check_fg("Primary amide", product_smiles)
                    or checker.check_fg("Secondary amide", product_smiles)
                    or checker.check_fg("Tertiary amide", product_smiles)
                )

                print(f"Product contains amide: {has_amide_in_product}")

                # Check for N-alkylation pattern which might not be detected as amide
                has_tertiary_amine_in_product = checker.check_fg("Tertiary amine", product_smiles)
                print(f"Product contains tertiary amine: {has_tertiary_amine_in_product}")

                # Check if this is an amide coupling reaction
                is_amide_coupling = (
                    checker.check_reaction(
                        "Acylation of Nitrogen Nucleophiles by Carboxylic Acids", rsmi
                    )
                    or checker.check_reaction(
                        "Acyl chloride with primary amine to amide (Schotten-Baumann)", rsmi
                    )
                    or checker.check_reaction("Acyl chloride with secondary amine to amide", rsmi)
                    or checker.check_reaction("Carboxylic acid with primary amine to amide", rsmi)
                    or checker.check_reaction("Ester with primary amine to amide", rsmi)
                    or checker.check_reaction("Ester with secondary amine to amide", rsmi)
                    or checker.check_reaction("Acyl chloride with ammonia to amide", rsmi)
                    or checker.check_reaction("Ester with ammonia to amide", rsmi)
                    or checker.check_reaction("Schotten-Baumann_amide", rsmi)
                    or checker.check_reaction(
                        "N-alkylation of primary amines with alkyl halides", rsmi
                    )
                    or checker.check_reaction(
                        "N-alkylation of secondary amines with alkyl halides", rsmi
                    )
                )

                print(f"Is known amide coupling reaction: {is_amide_coupling}")

                if is_amide_coupling:
                    print("Confirmed amide coupling reaction")
                    has_late_stage_amide = True
                elif has_amide_in_product:
                    # Check if amide is newly formed
                    # Count amides in reactants
                    amide_count_in_reactants = sum(
                        1
                        for r in reactants_smiles
                        if (
                            checker.check_fg("Primary amide", r)
                            or checker.check_fg("Secondary amide", r)
                            or checker.check_fg("Tertiary amide", r)
                        )
                    )

                    print(f"Amide count in reactants: {amide_count_in_reactants}")

                    # Check for carboxylic acid/acyl halide and amine in reactants
                    has_carboxylic_acid = any(
                        checker.check_fg("Carboxylic acid", r) for r in reactants_smiles
                    )
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

                    print(f"Has carboxylic acid: {has_carboxylic_acid}")
                    print(f"Has acyl halide: {has_acyl_halide}")
                    print(f"Has ester: {has_ester}")
                    print(f"Has primary amine: {has_primary_amine}")
                    print(f"Has secondary amine: {has_secondary_amine}")

                    # If product has amide and reactants have appropriate functional groups
                    # and not all reactants have amides, then amide is likely formed in this step
                    if amide_count_in_reactants < len(reactants_smiles) and (
                        (has_carboxylic_acid or has_acyl_halide or has_ester)
                        and (has_primary_amine or has_secondary_amine)
                    ):
                        print("Amide appears to be newly formed in final step")
                        has_late_stage_amide = True
                else:
                    # Check for N-alkylation pattern (carboxylic acid + amine → tertiary amine)
                    # This handles cases where the product is not detected as an amide but is chemically an amide coupling
                    has_carboxylic_acid = any(
                        checker.check_fg("Carboxylic acid", r) for r in reactants_smiles
                    )
                    has_primary_amine = any(
                        checker.check_fg("Primary amine", r) for r in reactants_smiles
                    )
                    has_secondary_amine = any(
                        checker.check_fg("Secondary amine", r) for r in reactants_smiles
                    )

                    print(f"Checking N-alkylation pattern:")
                    print(f"Has carboxylic acid: {has_carboxylic_acid}")
                    print(f"Has primary amine: {has_primary_amine}")
                    print(f"Has secondary amine: {has_secondary_amine}")
                    print(f"Has tertiary amine in product: {has_tertiary_amine_in_product}")

                    if (
                        has_tertiary_amine_in_product
                        and has_carboxylic_acid
                        and (has_primary_amine or has_secondary_amine)
                    ):
                        print("Detected N-alkylation pattern that represents amide coupling")
                        has_late_stage_amide = True

        # Continue traversing
        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    # Start traversal
    dfs_traverse(route)

    print(f"Final result - Detected late-stage amide coupling strategy: {has_late_stage_amide}")

    return has_late_stage_amide
