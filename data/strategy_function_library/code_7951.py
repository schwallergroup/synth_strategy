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

def main(route) -> Tuple[bool, Dict]:
    """
    This function detects if the synthesis route preserves both fluorine atoms
    and a cyano group throughout multiple steps.
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

    # Track paths with preserved functional groups
    preservation_paths = []
    result = False

    def dfs_traverse(node, current_path=None, depth=0):
        nonlocal findings_json
        if current_path is None:
            current_path = []

        # For molecule nodes, check for fluorine and cyano groups
        if node["type"] == "mol" and node.get("smiles"):
            mol_smiles = node["smiles"]
            has_fluorine = checker.check_fg("Trifluoro group", mol_smiles)
            if has_fluorine and "Trifluoro group" not in findings_json["atomic_checks"]["functional_groups"]:
                findings_json["atomic_checks"]["functional_groups"].append("Trifluoro group")

            has_cyano = checker.check_fg("Nitrile", mol_smiles)
            if has_cyano and "Nitrile" not in findings_json["atomic_checks"]["functional_groups"]:
                findings_json["atomic_checks"]["functional_groups"].append("Nitrile")

            # Create or update the current path information
            current_info = {
                "depth": depth,
                "has_fluorine": has_fluorine,
                "has_cyano": has_cyano,
                "smiles": mol_smiles,
            }

            # Add to the current path
            new_path = current_path + [current_info]

            # If this is a leaf node (starting material), save the path
            if node.get("in_stock", False) or not node.get("children"):
                if len(new_path) > 1:  # Only save paths with multiple nodes
                    preservation_paths.append(new_path)

            # Continue traversal with updated path
            for child in node.get("children", []):
                # Depth increases when going from chemical to reaction
                dfs_traverse(child, new_path, depth + 1)

        # For reaction nodes, just continue traversal
        elif node["type"] == "reaction":
            for child in node.get("children", []):
                # Depth remains the same when going from reaction to chemical
                dfs_traverse(child, current_path, depth)

    # Start traversal from the root
    dfs_traverse(route)

    # Check each path for preservation of both groups for at least 3 consecutive steps
    for path in preservation_paths:
        # Sort by depth to ensure correct order (target to starting material)
        path.sort(key=lambda x: x["depth"])

        # Count consecutive steps with both groups preserved
        consecutive_count = 0
        max_consecutive = 0

        for node_info in path:
            if node_info["has_fluorine"] and node_info["has_cyano"]:
                consecutive_count += 1
                max_consecutive = max(max_consecutive, consecutive_count)
            else:
                consecutive_count = 0

        if max_consecutive >= 3:
            print(f"Found preservation path with {max_consecutive} consecutive steps")
            result = True
            # Add the structural constraint to findings_json
            if {"type": "count", "details": {"target": "consecutive_molecule_steps_with_co-occurrence(Trifluoro group, Nitrile)", "operator": ">=", "value": 3}} not in findings_json["structural_constraints"]:
                findings_json["structural_constraints"].append({
                    "type": "count",
                    "details": {
                        "target": "consecutive_molecule_steps_with_co-occurrence(Trifluoro group, Nitrile)",
                        "operator": ">=",
                        "value": 3
                    }
                })
            # Since we found one path, we can break and return True
            break

    return result, findings_json