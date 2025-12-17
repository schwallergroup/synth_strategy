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
    Detects if the final product contains a cyclopropane ring that was also present in one of the starting materials. This identifies a building-block strategy where a cyclopropyl fragment is carried through the synthesis.
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

    max_depth = 0
    cyclopropyl_incorporation = False

    # First pass to determine the maximum depth of the tree
    def get_max_depth(node, current_depth=0):
        nonlocal max_depth
        max_depth = max(max_depth, current_depth)

        for child in node.get("children", []):
            get_max_depth(child, current_depth + 1)

    get_max_depth(route)

    # Define early stage threshold (e.g., first third of synthesis depth)
    early_stage_threshold = max_depth // 3

    # Second pass to check for cyclopropyl incorporation
    def dfs_traverse(node, depth=0):
        nonlocal cyclopropyl_incorporation, findings_json

        # Check if this is a molecule node
        if node["type"] == "mol":
            # Check if this molecule contains a cyclopropyl group
            if checker.check_ring("cyclopropane", node["smiles"]):
                if "cyclopropane" not in findings_json["atomic_checks"]["ring_systems"]:
                    findings_json["atomic_checks"]["ring_systems"].append("cyclopropane")

                # If this is the final product (depth 0) and it contains cyclopropyl
                if depth == 0:
                    # Check if any starting material had cyclopropyl
                    has_cyclopropyl_starting_material = False

                    def check_starting_materials(n):
                        nonlocal has_cyclopropyl_starting_material, findings_json
                        if n["type"] == "mol" and n.get("in_stock", False):
                            if checker.check_ring("cyclopropane", n["smiles"]):
                                if "cyclopropane" not in findings_json["atomic_checks"]["ring_systems"]:
                                    findings_json["atomic_checks"]["ring_systems"].append("cyclopropane")
                                has_cyclopropyl_starting_material = True
                        for c in n.get("children", []):
                            check_starting_materials(c)

                    check_starting_materials(route)

                    if has_cyclopropyl_starting_material:
                        cyclopropyl_incorporation = True
                        # Add the structural constraint if both conditions are met
                        findings_json["structural_constraints"].append(
                            {
                                "type": "co-occurrence",
                                "details": {
                                    "targets": [
                                        {
                                            "type": "positional",
                                            "details": {
                                                "target": "cyclopropane",
                                                "position": "final_product"
                                            }
                                        },
                                        {
                                            "type": "positional",
                                            "details": {
                                                "target": "cyclopropane",
                                                "position": "starting_material"
                                            }
                                        }
                                    ]
                                }
                            }
                        )

        # Continue traversal
        for child in node.get("children", []):
            # Determine the new depth based on the current node's type
            new_depth = depth
            if node["type"] != "reaction":  # Depth increases if not a reaction node
                new_depth = depth + 1
            dfs_traverse(child, new_depth)

    dfs_traverse(route)

    return cyclopropyl_incorporation, findings_json
