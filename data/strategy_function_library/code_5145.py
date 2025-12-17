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
    Detects a synthetic strategy where a trifluoromethyl group is preserved throughout the synthesis.
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

    # Track molecules with CF3 groups
    molecules_with_cf3 = []
    starting_materials_with_cf3 = []
    all_depths = set()
    result = True # Initialize the main boolean result

    def dfs_traverse(node, depth=0):
        nonlocal result # Declare result as nonlocal to modify it
        if node["type"] == "mol" and "smiles" in node:
            # Check if molecule has CF3 group
            has_cf3 = checker.check_fg("Trifluoro group", node["smiles"])

            if has_cf3:
                print(f"Found molecule with CF3 at depth {depth}: {node['smiles']}")
                if "Trifluoro group" not in findings_json["atomic_checks"]["functional_groups"]:
                    findings_json["atomic_checks"]["functional_groups"].append("Trifluoro group")
                molecules_with_cf3.append((depth, node["smiles"]))
                # Check if this is a starting material (either marked as in_stock or a leaf node)
                if node.get("in_stock", False) or (not node.get("children", [])):
                    print(f"Found starting material with CF3: {node['smiles']}")
                    starting_materials_with_cf3.append(node["smiles"])

        elif node["type"] == "reaction" and "metadata" in node and "mapped_reaction_smiles" in node["metadata"]:
            # Check if reaction preserves CF3 group
            rsmi = node["metadata"]["mapped_reaction_smiles"]
            reactants = rsmi.split(">")[0]
            products = rsmi.split(">")[-1]

            # If reactants have CF3, products should also have CF3
            reactants_has_cf3 = checker.check_fg("Trifluoro group", reactants)
            products_has_cf3 = checker.check_fg("Trifluoro group", products)

            if reactants_has_cf3 != products_has_cf3:
                print(f"CF3 group status changed in reaction: {rsmi}")
                result = False # Set result to False if CF3 status changes
                # Add the negation constraint if CF3 is modified
                if {"type": "negation", "details": {"target": "Trifluoro group modification in reaction"}} not in findings_json["structural_constraints"]:
                    findings_json["structural_constraints"].append({"type": "negation", "details": {"target": "Trifluoro group modification in reaction"}})
                return False

        # Determine the next depth based on the current node's type
        next_depth = depth
        if node["type"] != "reaction": # If current node is 'mol', depth increases for children
            next_depth = depth + 1

        # Traverse children
        for child in node.get("children", []):
            child_result = dfs_traverse(child, next_depth)
            if child_result is False:  # If any child traversal returns False, propagate it up
                result = False # Propagate the False result
                return False

        return True  # This branch of traversal is valid

    # Then traverse and check CF3 preservation
    if not dfs_traverse(route):
        result = False

    # Check if we found any molecules with CF3
    if not molecules_with_cf3:
        print("No molecules with CF3 groups found")
        result = False

    # Check if at least one starting material has CF3
    if not starting_materials_with_cf3:
        print("No starting materials with CF3 groups found")
        result = False
    else:
        # Add the count constraint if starting materials with CF3 are found
        if {"type": "count", "details": {"target": "starting_material_with_Trifluoro_group", "operator": ">", "value": 0}} not in findings_json["structural_constraints"]:
            findings_json["structural_constraints"].append({"type": "count", "details": {"target": "starting_material_with_Trifluoro_group", "operator": ">", "value": 0}})

    # Check if final product has CF3 (depth 0)
    final_product_has_cf3 = any(depth == 0 for depth, _ in molecules_with_cf3)
    if not final_product_has_cf3:
        print("Final product does not have CF3 group")
        result = False
    else:
        # Add the positional constraint if final product has CF3
        if {"type": "positional", "details": {"target": "Trifluoro group", "position": "last_stage"}} not in findings_json["structural_constraints"]:
            findings_json["structural_constraints"].append({"type": "positional", "details": {"target": "Trifluoro group", "position": "last_stage"}})

    if result:
        print("Trifluoromethyl group is preserved throughout synthesis")

    return result, findings_json
