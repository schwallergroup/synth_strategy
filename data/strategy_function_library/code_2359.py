from typing import Tuple, Dict, List
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


PERSISTENT_FGS_OF_INTEREST = ["Nitro group", "Allyl", "Tertiary amide"]

def main(route) -> Tuple[bool, Dict]:
    """
    This function detects if specific key functional groups are maintained throughout at least one complete synthetic path. It checks for the persistence of the following groups: Nitro group, Allyl, and Tertiary amide. A synthesis is flagged if at least two of these groups, present in the final product, are also present in all intermediate molecules back to the starting materials for any given path.
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

    # Define the functional groups to track
    target_functional_groups = PERSISTENT_FGS_OF_INTEREST

    # Track which functional groups are present in the target molecule
    target_fg_present = {fg: False for fg in target_functional_groups}

    # Process the target molecule first (depth 0)
    if route["type"] == "mol":
        mol_smiles = route["smiles"]
        for fg in target_functional_groups:
            if checker.check_fg(fg, mol_smiles):
                target_fg_present[fg] = True
                if fg not in findings_json["atomic_checks"]["functional_groups"]:
                    findings_json["atomic_checks"]["functional_groups"].append(fg)

    # Store complete paths where each functional group persists
    complete_paths = []

    # Track functional groups through each path
    def dfs_traverse(node, current_path, path_fg_present, depth=0):
        nonlocal findings_json
        # Make a copy of the current state to avoid modifying the parent's state
        current_path = current_path.copy()
        path_fg_present = copy.deepcopy(path_fg_present)

        # Add current node to path
        current_path.append(node)

        if node["type"] == "mol":
            mol_smiles = node["smiles"]

            # Check which functional groups are present in this molecule
            for fg in target_functional_groups:
                has_fg = checker.check_fg(fg, mol_smiles)
                if has_fg:
                    if fg not in findings_json["atomic_checks"]["functional_groups"]:
                        findings_json["atomic_checks"]["functional_groups"].append(fg)

                # Update the current path state - only track FGs that were in target
                if target_fg_present[fg]:
                    path_fg_present[fg] = has_fg

            # If this is a starting material (leaf node), we've completed a path
            if node.get("in_stock", False) or not node.get("children", []):
                # Store this complete path and its FG persistence
                complete_paths.append((current_path, path_fg_present))
                return

        # Determine the next depth based on the current node's type
        next_depth = depth
        if node["type"] != "reaction": # Depth increases if not a reaction node (i.e., chemical node)
            next_depth = depth + 1

        # Continue traversing children
        for child in node.get("children", []):
            dfs_traverse(child, current_path, path_fg_present, next_depth)

    # Start traversal from the root
    dfs_traverse(route, [], target_fg_present)

    # Count functional groups that persist in at least one complete path
    persistent_fgs = set()
    for path, fg_status in complete_paths:
        for fg in target_functional_groups:
            if target_fg_present[fg] and fg_status[fg]:
                persistent_fgs.add(fg)

    persistent_fg_count = len(persistent_fgs)

    result = persistent_fg_count >= 2

    # Add structural constraint if the condition is met
    if result:
        findings_json["structural_constraints"].append({
            "type": "count",
            "details": {
                "target": "persistent_functional_group",
                "operator": ">=",
                "value": 2,
                "condition": "A functional group from the candidate list is considered persistent if it is present in the final product and also in a starting material on at least one synthetic path."
            }
        })

    # Return True if at least two functional groups persist throughout
    return result, findings_json
