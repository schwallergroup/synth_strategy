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
from steerable_retro.utils.check import Check
from steerable_retro.utils import fuzzy_dict, check

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


ESTER_TO_AMIDE_REACTIONS = [
    "Aminolysis of esters",
    "Ester with primary amine to amide",
    "Ester with secondary amine to amide",
    "Ester with ammonia to amide",
]

def main(route) -> Tuple[bool, Dict]:
    """
    This function detects the use of specific ester-to-amide interconversion reactions during the later stages of a synthesis. It specifically checks for reactions defined in the `ESTER_TO_AMIDE_REACTIONS` list, including 'Aminolysis of esters', 'Ester with primary amine to amide', 'Ester with secondary amine to amide', and 'Ester with ammonia to amide'. A reaction is considered 'late-stage' if it occurs in the latter half of the synthetic route (i.e., when the step depth is less than or equal to half of the maximum depth).
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

    # First pass to determine the maximum depth
    max_depth = 0

    def calculate_max_depth(node, depth=0):
        nonlocal max_depth
        max_depth = max(max_depth, depth)
        for child in node.get("children", []):
            # Depth increases only when going from chemical to reaction
            # or from reaction to chemical (old logic for max_depth calculation)
            calculate_max_depth(child, depth + 1)

    calculate_max_depth(route)

    # Second pass to detect ester to amide conversions
    ester_to_amide_detected = False

    def dfs_traverse(node, depth=0):
        nonlocal ester_to_amide_detected, findings_json

        if node["type"] == "reaction":
            rsmi = node["metadata"].get("rsmi", "")
            if not rsmi:
                return

            # Check if the reaction is a known, specific ester-to-amide conversion
            is_ester_to_amide_reaction = False
            for reaction_name in ESTER_TO_AMIDE_REACTIONS:
                if checker.check_reaction(reaction_name, rsmi):
                    is_ester_to_amide_reaction = True
                    findings_json["atomic_checks"]["named_reactions"].append(reaction_name)
                    break

            if is_ester_to_amide_reaction:
                # Consider it late-stage if it's in the first half of the synthesis steps
                # (depth=1 is the final step)
                if depth <= max_depth / 2:
                    ester_to_amide_detected = True
                    # Add the structural constraint if the condition is met
                    findings_json["structural_constraints"].append({
                        "type": "positional",
                        "details": {
                            "targets": [
                                "Aminolysis of esters",
                                "Ester with primary amine to amide",
                                "Ester with secondary amine to amide",
                                "Ester with ammonia to amide"
                            ],
                            "position": "latter_half"
                        }
                    })

        # Traverse children
        for child in node.get("children", []):
            # New logic for depth calculation:
            # Depth increases only when traversing from a chemical node to a reaction node.
            # Depth remains the same when traversing from a reaction node to a chemical node.
            new_depth = depth
            if node["type"] != "reaction": # This means current node is 'chemical'
                new_depth = depth + 1
            
            dfs_traverse(child, new_depth)

    # Start traversal from root
    dfs_traverse(route)

    return ester_to_amide_detected, findings_json
