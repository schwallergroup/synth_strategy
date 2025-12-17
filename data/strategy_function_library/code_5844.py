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


REDUCTIVE_AMINATION_TYPES = [
    "Reductive amination with aldehyde",
    "Reductive amination with ketone",
    "Reductive amination with alcohol",
]

def main(route) -> Tuple[bool, Dict]:
    """
    Checks for specific reductive amination reactions in the final synthetic step (depth=1),
    including 'Reductive amination with aldehyde', 'Reductive amination with ketone', and 'Reductive amination with alcohol'.
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

    reductive_amination_detected = False

    def dfs_traverse(node, depth):
        nonlocal reductive_amination_detected, findings_json
        if reductive_amination_detected:  # Stop traversal if already detected
            return

        if node["type"] == "reaction" and depth == 1:
            rsmi = node["metadata"].get("rsmi", "")
            if rsmi:
                for rt in REDUCTIVE_AMINATION_TYPES:
                    if checker.check_reaction(rt, rsmi):
                        reductive_amination_detected = True
                        findings_json["atomic_checks"]["named_reactions"].append(rt)
                        # Since we found one, we can break from this inner loop
                        # but continue to allow other types to be added if they also match
                        # if the goal was to find *any* and stop, the outer if would handle it.
                        # For now, we append all found types.

        # Traverse children
        for child in node.get("children", []):
            # New logic for depth calculation
            new_depth = depth
            if node["type"] != "reaction": # If current node is chemical, depth increases
                new_depth = depth + 1
            # If current node is reaction, depth remains the same
            
            dfs_traverse(child, new_depth)

    # Start traversal from the root (product molecule at depth 0)
    dfs_traverse(route, 0)

    if reductive_amination_detected:
        # Add the structural constraint if any reductive amination was detected at depth 1
        findings_json["structural_constraints"].append({
            "type": "positional",
            "details": {
                "target": [
                    "Reductive amination with aldehyde",
                    "Reductive amination with ketone",
                    "Reductive amination with alcohol"
                ],
                "position": "final_step"
            }
        })

    return reductive_amination_detected, findings_json
