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


SUZUKI_REACTION_TYPES = [
    "Suzuki coupling with boronic acids",
    "Suzuki coupling with boronic acids OTf",
    "Suzuki coupling with sulfonic esters",
    "Suzuki coupling with boronic esters OTf",
    "Suzuki coupling with boronic esters",
    "{Suzuki}",
]

def main(route) -> Tuple[bool, Dict]:
    """
    Detects a late-stage Suzuki coupling reaction in the final synthetic step. This is
    identified by checking for specific named reaction types, including those defined
    in SUZUKI_REACTION_TYPES.
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

    final_step_is_suzuki = False

    def dfs_traverse(node, depth=0):
        nonlocal final_step_is_suzuki, findings_json

        # The final step is at depth=1. We only check this step.
        if node["type"] == "reaction" and depth == 1:
            try:
                rsmi = node["metadata"]["mapped_reaction_smiles"]
                for reaction_type in SUZUKI_REACTION_TYPES:
                    if checker.check_reaction(reaction_type, rsmi):
                        final_step_is_suzuki = True
                        findings_json["atomic_checks"]["named_reactions"].append(reaction_type)
                        # Add structural constraint if this is the final step and a Suzuki reaction is found
                        findings_json["structural_constraints"].append({
                            "type": "positional",
                            "details": {
                                "target": "Suzuki coupling",
                                "position": "last_stage"
                            }
                        })
                        break
            except Exception as e:
                print(f"Error processing reaction node: {e}")

        # Continue traversing the entire tree as per the original structure
        for child in node.get("children", []):
            # New logic for depth calculation
            new_depth = depth
            if node["type"] != "reaction": # If current node is chemical, depth increases
                new_depth = depth + 1
            # If current node is reaction, depth remains the same
            dfs_traverse(child, new_depth)

    # Start traversal
    dfs_traverse(route)
    return final_step_is_suzuki, findings_json
