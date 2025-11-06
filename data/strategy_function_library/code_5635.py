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


def main(route) -> Tuple[bool, Dict]:
    """
    Detects if the synthesis involves building a structure containing a trifluoromethyl group
    that is present from early stages and preserved throughout the synthesis.
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

    # Track CF3 groups at different stages
    cf3_by_depth = {}
    final_product_cf3 = False
    max_depth = 0

    def dfs_traverse(node, depth=0):
        nonlocal cf3_by_depth, final_product_cf3, max_depth, findings_json

        max_depth = max(max_depth, depth)

        if node["type"] == "mol":
            mol_smiles = node["smiles"]
            has_cf3 = checker.check_fg("Trifluoro group", mol_smiles)

            if has_cf3:
                if "Trifluoro group" not in findings_json["atomic_checks"]["functional_groups"]:
                    findings_json["atomic_checks"]["functional_groups"].append("Trifluoro group")
                if depth not in cf3_by_depth:
                    cf3_by_depth[depth] = set()
                cf3_by_depth[depth].add(mol_smiles)

            if depth == 0:
                if has_cf3:
                    final_product_cf3 = True

        # Traverse children
        for child in node.get("children", []):
            # New logic for depth calculation
            if node["type"] == "reaction":
                dfs_traverse(child, depth)
            else:
                dfs_traverse(child, depth + 1)

    # Start traversal from the root
    dfs_traverse(route)

    result = False

    if max_depth == 0:
        return result, findings_json

    # Define early stage as the deeper 70% of the synthesis
    early_stage_threshold = int(max_depth * 0.7)

    # Check if CF3 is present in early stages (deeper parts of the tree)
    early_cf3_present = any(depth >= early_stage_threshold for depth in cf3_by_depth.keys())
    if early_cf3_present:
        findings_json["structural_constraints"].append({
            "type": "positional",
            "details": {
                "target": "Trifluoro group",
                "position": "early_stage",
                "definition": "present in the deepest 70% of synthesis stages"
            }
        })

    # Check if CF3 is present throughout the synthesis
    depth_coverage = len(cf3_by_depth) / (max_depth + 1)
    cf3_throughout = depth_coverage >= 0.4
    if cf3_throughout:
        findings_json["structural_constraints"].append({
            "type": "count",
            "details": {
                "target": "stages_with_Trifluoro_group",
                "operator": ">=",
                "value": "40% of total stages"
            }
        })

    # Return true if CF3 is present from early stages and preserved to the final product
    result = early_cf3_present and final_product_cf3 and cf3_throughout

    if final_product_cf3:
        findings_json["structural_constraints"].append({
            "type": "positional",
            "details": {
                "target": "Trifluoro group",
                "position": "last_stage"
            }
        })

    return result, findings_json
