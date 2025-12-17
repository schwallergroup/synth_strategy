from typing import Tuple, Dict, List
import copy
from rdkit.Chem import AllChem, rdFMCS
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
from synth_strategy.utils.check import Check
from synth_strategy.utils import fuzzy_dict, check

from pathlib import Path
root_data = Path(__file__).parent.parent

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

def main(route) -> Tuple[bool, Dict]:
    """
    Detects a strategy involving late-stage functional group transformation,
    specifically looking for nitrile to carboxylic acid conversion in the final step.
    """
    findings_template = {
        "atomic_checks": {
            "named_reactions": [],
            "ring_systems": [],
            "functional_groups": []
        },
        "structural_constraints": []
    }
    findings_json = copy.deepcopy(findings_template)

    # Initialize tracking variables
    nitrile_hydrolysis_reactions = []
    all_reactions = []
    result = False

    def dfs_traverse(node, depth=0):
        nonlocal nitrile_hydrolysis_reactions, all_reactions, findings_json

        if node["type"] == "reaction":
            if "metadata" in node and "mapped_reaction_smiles" in node["metadata"]:
                rsmi = node["metadata"]["mapped_reaction_smiles"]
                reactants_part = rsmi.split(">")[0]
                product_part = rsmi.split(">")[-1]

                # Record all reactions with their depths
                all_reactions.append((depth, rsmi))

                # Check for nitrile hydrolysis using the checker functions
                has_nitrile = checker.check_fg("Nitrile", reactants_part)
                if has_nitrile:
                    if "Nitrile" not in findings_json["atomic_checks"]["functional_groups"]:
                        findings_json["atomic_checks"]["functional_groups"].append("Nitrile")

                has_carboxylic_acid = checker.check_fg("Carboxylic acid", product_part)
                if has_carboxylic_acid:
                    if "Carboxylic acid" not in findings_json["atomic_checks"]["functional_groups"]:
                        findings_json["atomic_checks"]["functional_groups"].append("Carboxylic acid")

                # Check if this is specifically a nitrile oxidation reaction
                is_nitrile_oxidation = checker.check_reaction(
                    "Oxidation of nitrile to carboxylic acid", rsmi
                )
                if is_nitrile_oxidation:
                    if "Oxidation of nitrile to carboxylic acid" not in findings_json["atomic_checks"]["named_reactions"]:
                        findings_json["atomic_checks"]["named_reactions"].append("Oxidation of nitrile to carboxylic acid")

                # If we can't directly identify the reaction, check for the functional group transformation
                # Make sure we're not just detecting a reaction where carboxylic acid was already present
                if (
                    has_nitrile
                    and has_carboxylic_acid
                    and (
                        is_nitrile_oxidation
                        or not checker.check_fg("Carboxylic acid", reactants_part)
                    )
                ):
                    print(f"Nitrile hydrolysis detected at depth {depth}")
                    nitrile_hydrolysis_reactions.append((depth, rsmi))

        # Continue traversing
        for child in node.get("children", []):
            # New logic for depth calculation
            new_depth = depth
            if node["type"] != "reaction": # If current node is chemical, depth increases
                new_depth = depth + 1
            # If current node is reaction, depth remains the same
            dfs_traverse(child, new_depth)

    # Start traversal
    dfs_traverse(route)

    # Sort reactions by depth (ascending)
    nitrile_hydrolysis_reactions.sort(key=lambda x: x[0])
    all_reactions.sort(key=lambda x: x[0])

    # Check if nitrile hydrolysis occurs in the final step
    final_step_nitrile_hydrolysis = False
    if nitrile_hydrolysis_reactions:
        min_depth_nitrile = nitrile_hydrolysis_reactions[0][0]
        print(f"Minimum depth of nitrile hydrolysis: {min_depth_nitrile}")

        # Find the minimum depth of any reaction
        min_depth_any = all_reactions[0][0] if all_reactions else float("inf")

        # Check if the nitrile hydrolysis is at the minimum depth (final step)
        # or at depth 1 which could be the final step in some routes
        if min_depth_nitrile == min_depth_any or min_depth_nitrile <= 1:
            # Count how many reactions are at this minimum depth
            reactions_at_min_depth = sum(1 for d, _ in all_reactions if d == min_depth_nitrile)

            # If there's only one reaction at this depth or all reactions at this depth are nitrile hydrolysis
            nitrile_reactions_at_min_depth = sum(
                1 for d, _ in nitrile_hydrolysis_reactions if d == min_depth_nitrile
            )

            if reactions_at_min_depth == nitrile_reactions_at_min_depth:
                final_step_nitrile_hydrolysis = True
                result = True
                # Add structural constraints if the condition is met
                findings_json["structural_constraints"].append({
                    "type": "positional",
                    "details": {
                        "target": "nitrile_to_carboxylic_acid_conversion",
                        "position": "last_stage"
                    }
                })
                findings_json["structural_constraints"].append({
                    "type": "negation",
                    "details": {
                        "target": "any_reaction_other_than_nitrile_conversion_in_last_stage"
                    }
                })

    return result, findings_json