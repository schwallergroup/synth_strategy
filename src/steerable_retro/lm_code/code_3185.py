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
    as the final step of the synthesis.
    """
    late_stage_amide_formation = False

    def dfs_traverse(node, depth=0):
        nonlocal late_stage_amide_formation

        print(f"Traversing node at depth {depth}, type: {node['type']}")

        if node["type"] == "reaction" and depth <= 1:  # Check at depths 0 and 1 (final steps)
            try:
                rsmi = node["metadata"]["rsmi"]
                print(f"Analyzing reaction SMILES: {rsmi}")

                reactants_smiles = rsmi.split(">")[0].split(".")
                product_smiles = rsmi.split(">")[-1]

                # Check if the product contains an amide group
                has_amide_product = (
                    checker.check_fg("Primary amide", product_smiles)
                    or checker.check_fg("Secondary amide", product_smiles)
                    or checker.check_fg("Tertiary amide", product_smiles)
                )

                print(f"Product contains amide: {has_amide_product}")

                if has_amide_product:
                    # Check for specific amide formation reactions
                    amide_reactions = [
                        "Acylation of Nitrogen Nucleophiles by Acyl/Thioacyl/Carbamoyl Halides and Analogs_N",
                        "Acylation of Nitrogen Nucleophiles by Acyl/Thioacyl/Carbamoyl Halides and Analogs_OS",
                        "Acyl chloride with ammonia to amide",
                        "Acyl chloride with primary amine to amide (Schotten-Baumann)",
                        "Acyl chloride with primary amine to imide",
                        "Acyl chloride with secondary amine to amide",
                        "Carboxylic acid with primary amine to amide",
                        "Ester with ammonia to amide",
                        "Ester with primary amine to amide",
                        "Ester with secondary amine to amide",
                        "Schotten-Baumann_amide",
                        "Acylation of primary amines",
                        "Acylation of secondary amines",
                        "Acylation of secondary amines with anhydrides",
                        "Carboxylic acid to amide conversion",
                        "Nitrile to amide",
                    ]

                    for reaction_type in amide_reactions:
                        if checker.check_reaction(reaction_type, rsmi):
                            print(
                                f"Detected late-stage amide formation as final step: {reaction_type}"
                            )
                            late_stage_amide_formation = True
                            return

                    # If no specific reaction type matched, check for reactant-product patterns
                    # First check if reactants don't already have the amide
                    reactants_have_amide = any(
                        checker.check_fg("Primary amide", r)
                        or checker.check_fg("Secondary amide", r)
                        or checker.check_fg("Tertiary amide", r)
                        for r in reactants_smiles
                    )

                    if not reactants_have_amide:
                        has_acyl_halide = any(
                            checker.check_fg("Acyl halide", r) for r in reactants_smiles
                        )
                        has_carboxylic_acid = any(
                            checker.check_fg("Carboxylic acid", r) for r in reactants_smiles
                        )
                        has_ester = any(checker.check_fg("Ester", r) for r in reactants_smiles)
                        has_anhydride = any(
                            checker.check_fg("Anhydride", r) for r in reactants_smiles
                        )
                        has_nitrile = any(checker.check_fg("Nitrile", r) for r in reactants_smiles)

                        has_primary_amine = any(
                            checker.check_fg("Primary amine", r) for r in reactants_smiles
                        )
                        has_secondary_amine = any(
                            checker.check_fg("Secondary amine", r) for r in reactants_smiles
                        )
                        has_amine = has_primary_amine or has_secondary_amine

                        if (
                            has_acyl_halide
                            or has_carboxylic_acid
                            or has_ester
                            or has_anhydride
                            or has_nitrile
                        ) and has_amine:
                            print(
                                "Detected late-stage amide formation as final step (reactant pattern)"
                            )
                            late_stage_amide_formation = True
                            return
            except Exception as e:
                print(f"Error analyzing reaction: {e}")

        # Traverse children
        for child in node.get("children", []):
            if (
                not late_stage_amide_formation
            ):  # Stop traversal if we already found what we're looking for
                dfs_traverse(child, depth + 1)

    # Start traversal from the root
    dfs_traverse(route)

    return late_stage_amide_formation
