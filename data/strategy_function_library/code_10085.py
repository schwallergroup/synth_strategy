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
rng_smiles_args = {
    "file_path": f"{root_data}/patterns/chemical_rings_smiles.json",
    "value_field": "smiles",
    "key_field": "name",
}
functional_groups = fuzzy_dict.FuzzyDict.from_json(**fg_args)
reaction_classes = fuzzy_dict.FuzzyDict.from_json(**reaction_class_args)
ring_smiles = fuzzy_dict.FuzzyDict.from_json(**rng_smiles_args)

checker = check.Check(
    fg_dict=functional_groups, reaction_dict=reaction_classes, ring_dict=ring_smiles
)


def main(route) -> Tuple[bool, Dict]:
    """
    Detects if an amide functionality is maintained throughout the synthesis.
    Returns True if all non-starting molecules in the main synthetic pathway contain an amide group.
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

    # Track if all molecules in the main synthetic pathway contain an amide group
    all_have_amide = True

    def has_amide(smiles):
        """Check if a molecule has any type of amide group"""
        found_fgs = []
        if checker.check_fg("Primary amide", smiles):
            found_fgs.append("Primary amide")
        if checker.check_fg("Secondary amide", smiles):
            found_fgs.append("Secondary amide")
        if checker.check_fg("Tertiary amide", smiles):
            found_fgs.append("Tertiary amide")
        
        # Add found functional groups to findings_json, avoiding duplicates
        for fg in found_fgs:
            if fg not in findings_json["atomic_checks"]["functional_groups"]:
                findings_json["atomic_checks"]["functional_groups"].append(fg)

        return len(found_fgs) > 0

    def dfs_traverse(node, depth=0, is_main_path=True):
        nonlocal all_have_amide

        if node["type"] == "mol":
            # Only check non-starting materials in the main synthetic pathway
            if not node.get("in_stock", False) and is_main_path:
                smiles = node.get("smiles", "")
                if smiles:
                    if not has_amide(smiles):
                        print(f"Molecule at depth {depth} lacks amide functionality: {smiles}")
                        all_have_amide = False
                    else:
                        print(f"Molecule at depth {depth} has amide functionality: {smiles}")

            # Continue traversal
            for child in node.get("children", []):
                # If this is a molecule node, all children are reaction nodes
                # and should remain in the main path
                # Depth increases when going from chemical (mol) to reaction
                dfs_traverse(child, depth + 1, is_main_path)

        elif node["type"] == "reaction":
            # For reaction nodes, identify the main product
            # The main product is typically the first child (in retrosynthetic direction)
            children = node.get("children", [])

            if children:
                # The first child is the main product in retrosynthetic direction
                # Depth remains the same when going from reaction to chemical
                main_product = children[0]
                dfs_traverse(main_product, depth, True)

                # Other children are reagents/catalysts, not in the main path
                # Depth remains the same when going from reaction to chemical
                for child in children[1:]:
                    dfs_traverse(child, depth, False)

    # First check the target molecule (root of the tree)
    result = True
    if route["type"] == "mol" and not route.get("in_stock", False):
        smiles = route.get("smiles", "")
        if smiles and not has_amide(smiles):
            print(f"Target molecule lacks amide functionality: {smiles}")
            result = False
        else:
            print(f"Target molecule has amide functionality: {smiles}")

    # Then traverse the rest of the tree
    dfs_traverse(route)

    if all_have_amide and result:
        print("Detected maintenance of amide functionality throughout synthesis")
        # Add the structural constraint if the condition is met
        findings_json["structural_constraints"].append({
            "type": "count",
            "details": {
                "target": "intermediate_lacking_amide_group",
                "operator": "==",
                "value": 0
            }
        })
    else:
        print("Amide functionality is not maintained throughout synthesis")
        result = False # Ensure result is False if all_have_amide is False

    return result and all_have_amide, findings_json
