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
    This function detects a synthetic strategy where a thiophene core is preserved
    throughout the synthesis. The thiophene ring should be present in the final product
    and all intermediate products, but not necessarily in reagents or starting materials.
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

    # Track if the main synthetic pathway maintains a thiophene core
    thiophene_in_main_pathway = True
    target_has_thiophene = False

    def dfs_traverse(node, depth=0, is_main_product=True):
        nonlocal thiophene_in_main_pathway, findings_json

        # Only check molecules that are in the main product line
        if node["type"] == "mol" and node.get("smiles") and is_main_product:
            mol_smiles = node["smiles"]
            is_starting_material = node.get("in_stock", False)

            # Check if molecule contains thiophene, but only if it's not a starting material
            if not is_starting_material:
                has_thiophene = checker.check_ring("thiophene", mol_smiles)
                if has_thiophene:
                    if "thiophene" not in findings_json["atomic_checks"]["ring_systems"]:
                        findings_json["atomic_checks"]["ring_systems"].append("thiophene")
                else:
                    thiophene_in_main_pathway = False

        # Determine the next depth based on the current node's type
        next_depth = depth
        if node["type"] != "reaction":  # If current node is chemical, depth increases
            next_depth = depth + 1

        # Process children nodes
        if node["type"] == "reaction" and node.get("children"):
            # In retrosynthetic traversal, the first child is the product
            for i, child in enumerate(node.get("children", [])):
                # Only the first child is considered the main product in retrosynthetic direction
                dfs_traverse(child, next_depth, is_main_product=(i == 0))
        else:
            # Process children for molecule nodes (maintaining the is_main_product flag)
            for child in node.get("children", []):
                dfs_traverse(child, next_depth, is_main_product)

    result = False
    # Start traversal with the target molecule
    try:
        # Check if the target molecule contains thiophene
        if route["type"] == "mol" and route.get("smiles"):
            target_has_thiophene = checker.check_ring("thiophene", route["smiles"])
            if target_has_thiophene:
                if "thiophene" not in findings_json["atomic_checks"]["ring_systems"]:
                    findings_json["atomic_checks"]["ring_systems"].append("thiophene")
            else:
                result = False
                return result, findings_json

        # Traverse the synthesis route
        dfs_traverse(route)

        result = thiophene_in_main_pathway

        if result and target_has_thiophene:
            findings_json["structural_constraints"].append({
                "type": "co-occurrence",
                "details": {
                    "targets": [
                        "thiophene_in_target_molecule",
                        "thiophene_in_all_main_path_intermediates"
                    ]
                }
            })

        return result, findings_json
    except Exception:
        return False, findings_json
