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
    Detects synthesis routes where a specific fluorinated functional group
    (Trifluoro group or Triflate) is installed early and retained throughout the synthesis.
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

    # List of fluorinated functional groups to check
    fluorinated_groups = ["Trifluoro group", "Triflate"]

    # Track fluorinated groups and their paths
    paths_with_fluorinated_groups = []
    max_depth = 0

    def has_fluorinated_group(smiles):
        """Check if molecule contains any fluorinated group and return the group names"""
        if not smiles:
            return []

        found_groups = []

        # Check for explicitly fluorinated groups
        for fg in fluorinated_groups:
            if checker.check_fg(fg, smiles):
                found_groups.append(fg)
                if fg not in findings_json["atomic_checks"]["functional_groups"]:
                    findings_json["atomic_checks"]["functional_groups"].append(fg)

        return found_groups

    def dfs_traverse(node, depth=0, current_path=None, fluorinated_at_depth=None):
        nonlocal max_depth
        max_depth = max(max_depth, depth)

        if current_path is None:
            current_path = []

        if fluorinated_at_depth is None:
            fluorinated_at_depth = {}

        current_path.append((node, depth))

        if node["type"] == "mol":
            smiles = node["smiles"]
            if smiles:
                fluorinated_groups_found = has_fluorinated_group(smiles)
                fluorinated_at_depth[depth] = fluorinated_groups_found

                if not node.get("children", []) or node.get("in_stock", False):
                    if any(fluorinated_at_depth.values()):
                        paths_with_fluorinated_groups.append(
                            (current_path.copy(), fluorinated_at_depth.copy())
                        )

        for child in node.get("children", []):
            # Determine the new depth based on the current node's type
            new_depth = depth
            if node["type"] != "reaction": # Depth increases only when not traversing from a reaction node
                new_depth = depth + 1
            
            dfs_traverse(child, new_depth, current_path.copy(), fluorinated_at_depth.copy())

    dfs_traverse(route)

    result = False
    if not paths_with_fluorinated_groups:
        return result, findings_json

    valid_paths = []
    for path, fluorinated_depths in paths_with_fluorinated_groups:
        depths_with_fluorine = [d for d, groups in fluorinated_depths.items() if groups]
        if not depths_with_fluorine:
            continue

        introduction_depth = max(depths_with_fluorine)

        retained = True
        loss_detected = False
        for i in range(len(path) - 1):
            current_node, current_depth = path[i]
            next_node, next_depth = path[i + 1]

            if current_depth < introduction_depth:
                continue

            # Corrected logic: check for loss (present in reactant, absent in product)
            if (
                next_depth in fluorinated_depths
                and fluorinated_depths[next_depth]
                and current_depth in fluorinated_depths
                and not fluorinated_depths[current_depth]
            ):
                retained = False
                loss_detected = True
                break

        if retained and introduction_depth > 1:
            valid_paths.append((path, introduction_depth))
            # Record structural constraints if conditions are met
            if {"type": "positional", "details": {"target": "introduction_of_fluorinated_group", "position": "depth > 1"}} not in findings_json["structural_constraints"]:
                findings_json["structural_constraints"].append({"type": "positional", "details": {"target": "introduction_of_fluorinated_group", "position": "depth > 1"}})
            # If retained, then negation of loss is true
            if {"type": "negation", "details": {"target": "loss_of_fluorinated_group", "condition": "after_introduction"}} not in findings_json["structural_constraints"]:
                findings_json["structural_constraints"].append({"type": "negation", "details": {"target": "loss_of_fluorinated_group", "condition": "after_introduction"}})
        elif loss_detected:
            # If loss was detected, the negation constraint is not met
            # We don't add the negation constraint if it was not met
            pass

    result = len(valid_paths) > 0
    return result, findings_json
