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
    This function detects a synthetic strategy involving incorporation of a triazole ring
    as a building block.
    """
    print("Starting triazole incorporation strategy analysis")

    findings_template = {
        "atomic_checks": {
            "named_reactions": [],
            "ring_systems": [],
            "functional_groups": []
        },
        "structural_constraints": []
    }
    findings_json = copy.deepcopy(findings_template)

    # Track if triazole is in final product and in any starting material
    triazole_in_final = False
    triazole_as_building_block = False

    def dfs_traverse(node, depth=0):
        nonlocal triazole_in_final, triazole_as_building_block, findings_json

        if node["type"] == "mol":
            mol_smiles = node["smiles"]
            print(f"Checking molecule at depth {depth}: {mol_smiles[:30]}...")

            try:
                # Check if molecule contains a triazole ring using the checker function
                has_triazole = checker.check_ring("triazole", mol_smiles)

                if has_triazole:
                    findings_json["atomic_checks"]["ring_systems"].append("triazole")
                    if depth == 0:  # Final product
                        print(f"Triazole ring detected in final product: {mol_smiles[:30]}...")
                        triazole_in_final = True
                        findings_json["structural_constraints"].append({"type": "positional", "details": {"target": "triazole", "position": "last_stage"}})
                    elif node.get("in_stock", False) or not node.get(
                        "children", []
                    ):  # Starting material
                        print(
                            f"Triazole ring detected in starting material at depth {depth}: {mol_smiles[:30]}... (in_stock: {node.get('in_stock', False)})"
                        )
                        triazole_as_building_block = True
                        findings_json["structural_constraints"].append({"type": "positional", "details": {"target": "triazole", "position": "starting_material"}})
            except Exception as e:
                print(f"Error processing molecule SMILES for triazole detection: {str(e)}")

        # Determine the next depth based on the current node's type
        next_depth = depth
        if node["type"] != "reaction":  # If current node is 'mol' (chemical), depth increases
            next_depth = depth + 1

        # Traverse children
        for child in node.get("children", []):
            dfs_traverse(child, next_depth)

    # Start traversal from the root
    dfs_traverse(route)

    # Return True if triazole is in final product and was incorporated as a building block
    result = triazole_in_final and triazole_as_building_block
    print(f"Triazole in final product: {triazole_in_final}")
    print(f"Triazole as building block: {triazole_as_building_block}")
    print(f"Triazole incorporation strategy detected: {result}")

    if result:
        findings_json["structural_constraints"].append({"type": "co-occurrence", "details": {"targets": ["triazole_at_last_stage", "triazole_at_starting_material"], "description": "A triazole ring must be present in the final product (last stage) and also be present in a starting material, indicating it was used as a building block."}})

    return result, findings_json
