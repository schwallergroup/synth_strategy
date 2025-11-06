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


SUZUKI_REACTION_NAMES = [
    "Suzuki coupling with boronic acids",
    "Suzuki coupling with boronic acids OTf",
    "Suzuki coupling with boronic esters",
    "Suzuki coupling with boronic esters OTf",
    "Suzuki coupling with sulfonic esters",
    "Suzuki",  # Generic Suzuki check
]

def main(route) -> Tuple[bool, Dict]:
    """
    Detects if the synthesis uses a late-stage Suzuki coupling (depth 0 or 1)
    for final diversification. This check identifies several Suzuki reaction types
    by name, as defined in SUZUKI_REACTION_NAMES.
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

    suzuki_detected = False
    max_depth = 0

    def dfs_traverse(node, depth=0):
        nonlocal suzuki_detected, max_depth, findings_json

        max_depth = max(max_depth, depth)
        print(f"Traversing node at depth {depth}, type: {node['type']}")

        if node["type"] == "reaction" and "metadata" in node and "mapped_reaction_smiles" in node["metadata"]:
            rsmi = node["metadata"]["mapped_reaction_smiles"]
            print(f"Checking reaction at depth {depth}: {rsmi}")

            # Check if this is a Suzuki coupling using the checker function
            is_suzuki = False
            for name in SUZUKI_REACTION_NAMES:
                if checker.check_reaction(name, rsmi):
                    is_suzuki = True
                    if name not in findings_json["atomic_checks"]["named_reactions"]:
                        findings_json["atomic_checks"]["named_reactions"].append(name)
                    break

            if is_suzuki:
                print(f"Detected Suzuki coupling at depth {depth}")
                if depth <= 1:
                    print("This is a late-stage Suzuki coupling")
                    suzuki_detected = True
                    # Add structural constraint if late-stage Suzuki is detected
                    if {"type": "positional", "details": {"target": "Suzuki coupling", "position": "late_stage"}} not in findings_json["structural_constraints"]:
                        findings_json["structural_constraints"].append({
                            "type": "positional",
                            "details": {
                                "target": "Suzuki coupling",
                                "position": "late_stage"
                            }
                        })

        for child in node.get("children", []):
            # New logic for depth calculation
            new_depth = depth
            if node['type'] != 'reaction': # If current node is chemical, depth increases
                new_depth = depth + 1
            # If current node is reaction, depth remains the same
            dfs_traverse(child, new_depth)

    dfs_traverse(route)
    print(f"Final result: suzuki_detected={suzuki_detected}, max_depth={max_depth}")

    return suzuki_detected, findings_json
