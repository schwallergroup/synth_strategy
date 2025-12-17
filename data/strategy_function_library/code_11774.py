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

def main(route) -> Tuple[bool, Dict]:
    """
    Detects a synthesis strategy where a trifluoromethyl group is present in a starting material
    and preserved through to the final product.
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

    # Track if we found the pattern
    cf3_in_final_product = False
    cf3_in_starting_material = False

    def dfs_traverse(node, depth=0):
        nonlocal cf3_in_final_product, cf3_in_starting_material, findings_json

        if node["type"] == "mol":
            # Check if this is a molecule node
            mol_smiles = node["smiles"]

            # Check for CF3 group using the checker function
            has_cf3 = checker.check_fg("Trifluoro group", mol_smiles)

            if has_cf3:
                if "Trifluoro group" not in findings_json["atomic_checks"]["functional_groups"]:
                    findings_json["atomic_checks"]["functional_groups"].append("Trifluoro group")

            # If this is the root node (final product)
            if depth == 0:
                cf3_in_final_product = has_cf3
                if cf3_in_final_product:
                    findings_json["structural_constraints"].append({
                        "type": "positional",
                        "details": {
                            "target": "Trifluoro group",
                            "position": "last_stage"
                        }
                    })
                print(f"Final product has CF3: {has_cf3}, SMILES: {mol_smiles}")

            # If this is a leaf node (starting material)
            is_leaf = len(node.get("children", [])) == 0 or node.get("in_stock", False)
            if is_leaf and has_cf3:
                cf3_in_starting_material = True
                findings_json["structural_constraints"].append({
                    "type": "positional",
                    "details": {
                        "target": "Trifluoro group",
                        "position": "starting_material"
                    }
                })
                print(f"Found CF3 in starting material at depth {depth}, SMILES: {mol_smiles}")

        # Traverse children
        for child in node.get("children", []):
            # Determine the new depth based on the current node's type
            new_depth = depth
            if node["type"] != "reaction": # If current node is 'mol', depth increases
                new_depth = depth + 1
            # If current node is 'reaction', depth remains the same
            
            dfs_traverse(child, new_depth)

    # Start traversal
    dfs_traverse(route)

    print(f"CF3 in final product: {cf3_in_final_product}")
    print(f"CF3 in starting material: {cf3_in_starting_material}")

    # Return True if CF3 is both in a starting material and the final product
    result = cf3_in_final_product and cf3_in_starting_material
    return result, findings_json
