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


FORWARD_REDUCTION_REACTIONS = [
    "Reduction of aldehydes and ketones to alcohols",
    "Reduction of ester to primary alcohol",
    "Reduction of ketone to secondary alcohol",
    "Reduction of carboxylic acid to primary alcohol",
    "Reduction of nitrile to amine",
    "Reduction of primary amides to amines",
    "Reduction of secondary amides to amines",
    "Reduction of tertiary amides to amines",
]

FORWARD_OXIDATION_REACTIONS = [
    "Oxidation of aldehydes to carboxylic acids",
    "Oxidation or Dehydrogenation of Alcohols to Aldehydes and Ketones",
    "Oxidation of alcohol to carboxylic acid",
    "Oxidation of ketone to carboxylic acid",
    "Oxidation of nitrile to carboxylic acid",
    "Oxidation of amide to carboxylic acid",
    "Oxidation of alkene to carboxylic acid",
    "Oxidation of alkene to aldehyde",
    "Oxidation of alcohol and aldehyde to ester",
    "Oxidation of boronic acids",
    "Oxidation of boronic esters",
    "Oxidative esterification of primary alcohols",
]

def main(route) -> Tuple[bool, Dict]:
    """
    Detects if the route employs a sequence of redox manipulations
    (reduction followed by oxidation or vice versa).
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

    # Track redox reactions with their depths
    redox_reactions = []

    def dfs_traverse(node, depth=0):
        nonlocal redox_reactions, findings_json
        if node["type"] == "reaction":
            try:
                rsmi = node["metadata"].get("rsmi", "")
                if not rsmi:
                    return

                # Check for reduction reactions (in forward direction)
                for r in FORWARD_REDUCTION_REACTIONS:
                    if checker.check_reaction(r, rsmi):
                        redox_reactions.append(("reduction", depth))
                        if r not in findings_json["atomic_checks"]["named_reactions"]:
                            findings_json["atomic_checks"]["named_reactions"].append(r)
                        break # Assuming only one match per reaction node is sufficient

                # Check for oxidation reactions (in forward direction)
                for r in FORWARD_OXIDATION_REACTIONS:
                    if checker.check_reaction(r, rsmi):
                        redox_reactions.append(("oxidation", depth))
                        if r not in findings_json["atomic_checks"]["named_reactions"]:
                            findings_json["atomic_checks"]["named_reactions"].append(r)
                        break # Assuming only one match per reaction node is sufficient

            except Exception as e:
                print(f"Error processing reaction at depth {depth}: {e}")

        # Continue traversal
        for child in node.get("children", []):
            # New logic: depth increases only when traversing from chemical to reaction
            new_depth = depth + 1 if node["type"] != "reaction" else depth
            dfs_traverse(child, new_depth)

    # Start traversal
    dfs_traverse(route)

    # Check if there's a sequence of redox manipulations
    has_redox_sequence = False

    # Sort by depth to make checking adjacent steps easier
    redox_reactions.sort(key=lambda x: x[1])

    # Check for redox steps of different types with at most one step in between
    for i in range(len(redox_reactions) - 1):
        current_type, current_depth = redox_reactions[i]
        next_type, next_depth = redox_reactions[i + 1]

        if current_type != next_type and abs(current_depth - next_depth) <= 2:
            has_redox_sequence = True
            # Add the structural constraint to findings_json
            findings_json["structural_constraints"].append({
                "type": "sequence",
                "details": {
                    "description": "Finds a sequence where a reduction reaction is followed by an oxidation reaction (or vice versa). The two reactions must be consecutive in the ordered list of all redox reactions in the route, and their depth difference must be at most 2.",
                    "event_A_group": {
                        "group_name": "reduction",
                        "members": [
                            "Reduction of aldehydes and ketones to alcohols",
                            "Reduction of ester to primary alcohol",
                            "Reduction of ketone to secondary alcohol",
                            "Reduction of carboxylic acid to primary alcohol",
                            "Reduction of nitrile to amine",
                            "Reduction of primary amides to amines",
                            "Reduction of secondary amides to amines",
                            "Reduction of tertiary amides to amines"
                        ]
                    },
                    "event_B_group": {
                        "group_name": "oxidation",
                        "members": [
                            "Oxidation of aldehydes to carboxylic acids",
                            "Oxidation or Dehydrogenation of Alcohols to Aldehydes and Ketones",
                            "Oxidation of alcohol to carboxylic acid",
                            "Oxidation of ketone to carboxylic acid",
                            "Oxidation of nitrile to carboxylic acid",
                            "Oxidation of amide to carboxylic acid",
                            "Oxidation of alkene to carboxylic acid",
                            "Oxidation of alkene to aldehyde",
                            "Oxidation of alcohol and aldehyde to ester",
                            "Oxidation of boronic acids",
                            "Oxidation of boronic esters",
                            "Oxidative esterification of primary alcohols"
                        ]
                    },
                    "max_depth_difference": 2
                }
            })
            break

    return has_redox_sequence, findings_json
