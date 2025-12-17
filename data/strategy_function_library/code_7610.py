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
    This function detects if a protected amine (specifically Boc-protected)
    is maintained throughout the synthesis without deprotection/reprotection.
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

    # Track protected amines at each depth and deprotection reactions
    protected_amines_by_depth = {}
    deprotection_depths = set()
    result = False

    def dfs_traverse(node, depth=0):
        nonlocal result, findings_json
        if node["type"] == "mol" and "smiles" in node:
            # Check for Boc-protected amine in molecule nodes
            if checker.check_fg("Boc", node["smiles"]):
                protected_amines_by_depth[depth] = node["smiles"]
                findings_json["atomic_checks"]["functional_groups"].append("Boc")
                # print(f"Found Boc-protected amine at depth {depth}: {node['smiles']}")

        elif node["type"] == "reaction" and "metadata" in node and "mapped_reaction_smiles" in node["metadata"]:
            # Check for Boc deprotection reactions
            if checker.check_reaction("Boc amine deprotection", node["metadata"]["mapped_reaction_smiles"]):
                deprotection_depths.add(depth)
                findings_json["atomic_checks"]["named_reactions"].append("Boc amine deprotection")
                # print(f"Found Boc deprotection at depth {depth}")

        # Traverse children
        for child in node.get("children", []):
            # New logic: depth increases only when traversing from chemical to reaction node
            # Depth remains the same when traversing from reaction to chemical node
            new_depth = depth
            if node["type"] != "reaction": # This means current node is 'mol' or 'chemical'
                new_depth = depth + 1
            dfs_traverse(child, new_depth)

    # Start traversal
    dfs_traverse(route)

    # Check if protected amine is present at multiple depths (maintained throughout)
    if len(protected_amines_by_depth) >= 2:
        min_depth = min(protected_amines_by_depth.keys())
        max_depth = max(protected_amines_by_depth.keys())

        # If protected amine is present at both early and late stages
        if max_depth - min_depth >= 2:
            # Add structural constraint for Boc presence at multiple stages
            findings_json["structural_constraints"].append({
                "type": "count",
                "details": {
                    "target": "Boc",
                    "operator": ">=",
                    "value": 2,
                    "note": "Counts the number of synthesis stages where a Boc group is present."
                }
            })
            # Check if there are no deprotection reactions ANYWHERE in the synthesis
            if not deprotection_depths:
                # print("Strategy detected: Protected amine maintained throughout synthesis")
                result = True
                # Add structural constraint for no deprotection
                findings_json["structural_constraints"].append({
                    "type": "negation",
                    "details": {
                        "target": "Boc amine deprotection",
                        "note": "The synthesis must not contain any Boc amine deprotection steps."
                    }
                })
            # else:
                # print("Boc group was deprotected during synthesis")

    return result, findings_json
