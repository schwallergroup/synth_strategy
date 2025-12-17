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
    This function detects a synthetic strategy where a tetrazole group is maintained
    throughout the synthesis (present in early and final stages).
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

    # Track tetrazole presence at different depths
    tetrazole_info = []

    def dfs_traverse(node, depth=0):
        nonlocal tetrazole_info, findings_json

        if node["type"] == "mol":
            # Check for tetrazole in molecule using the checker function
            if checker.check_ring("tetrazole", node["smiles"]):
                # Get atom indices of tetrazole ring
                tetrazole_indices = checker.get_ring_atom_indices("tetrazole", node["smiles"])
                tetrazole_info.append((depth, node["smiles"], tetrazole_indices))
                if "tetrazole" not in findings_json["atomic_checks"]["ring_systems"]:
                    findings_json["atomic_checks"]["ring_systems"].append("tetrazole")

        # Traverse children
        for child in node.get("children", []):
            # New logic: depth increases only when going from chemical to reaction node
            # Depth remains same when going from reaction to chemical node
            new_depth = depth
            if node["type"] != "reaction": # This means current node is 'mol' or 'chemical'
                new_depth = depth + 1
            
            dfs_traverse(child, new_depth)

    # Start traversal
    dfs_traverse(route)

    # No tetrazoles found
    if not tetrazole_info:
        return False, findings_json

    # Sort by depth to find early and late stages
    tetrazole_info.sort(key=lambda x: x[0])

    # Check if tetrazole is present in final product (depth 0)
    has_late = tetrazole_info[0][0] == 0
    if has_late:
        findings_json["structural_constraints"].append({
            "type": "positional",
            "details": {
                "target": "tetrazole",
                "position": "last_stage"
            }
        })

    # Find maximum depth (early stage)
    max_depth = tetrazole_info[-1][0]
    has_early = max_depth >= 2  # Early stage defined as depth >= 2
    if has_early:
        findings_json["structural_constraints"].append({
            "type": "positional",
            "details": {
                "target": "tetrazole",
                "position": "early_stage"
            }
        })

    # Check if tetrazole is maintained throughout the synthesis
    # This requires tetrazole to be present at both early and late stages
    result = has_early and has_late

    return result, findings_json