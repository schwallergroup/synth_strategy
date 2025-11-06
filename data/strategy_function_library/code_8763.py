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


HALOGEN_FGS = [
    "Primary halide",
    "Secondary halide",
    "Tertiary halide",
    "Aromatic halide",
    "Alkenyl halide",
    "Haloalkyne",
]

def main(route) -> Tuple[bool, Dict]:
    """
    This function detects if a halogen substituent is maintained throughout the synthesis.
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

    # Track halogen atoms through the synthesis
    halogen_maintained = False

    # Get the final product (root node)
    final_product_smiles = route["smiles"]
    final_product = Chem.MolFromSmiles(final_product_smiles)

    # Check if final product has any halogen
    has_halogen = False
    for halogen_type in HALOGEN_FGS:
        if checker.check_fg(halogen_type, final_product_smiles):
            has_halogen = True
            findings_json["atomic_checks"]["functional_groups"].append(halogen_type)
            print(f"Final product contains {halogen_type}")
            break

    if not has_halogen:
        print("Final product does not contain any halogen")
        return False, findings_json
    else:
        # If final product has halogen, add the positional constraint
        findings_json["structural_constraints"].append({
            "type": "positional",
            "details": {
                "target": [
                    "Primary halide",
                    "Secondary halide",
                    "Tertiary halide",
                    "Aromatic halide",
                    "Alkenyl halide",
                    "Haloalkyne"
                ],
                "position": "last_stage"
            }
        })

    # Track if halogen is maintained throughout synthesis
    def dfs_traverse(node, depth=0):
        nonlocal halogen_maintained, findings_json

        # Process molecule nodes
        if node["type"] == "mol":
            mol_smiles = node["smiles"]

            # Skip checking starting materials (in_stock)
            if node.get("in_stock", False):
                return

            # Check if this intermediate has a halogen
            has_halogen_intermediate = False
            for halogen_type in HALOGEN_FGS:
                if checker.check_fg(halogen_type, mol_smiles):
                    has_halogen_intermediate = True
                    findings_json["atomic_checks"]["functional_groups"].append(halogen_type)
                    print(f"Intermediate at depth {depth} contains {halogen_type}")
                    break

            if not has_halogen_intermediate and depth > 0:  # Skip final product (depth 0)
                print(f"Halogen not maintained in intermediate at depth {depth}")
                halogen_maintained = False
                # Add negation constraint if an intermediate lacks halogen
                findings_json["structural_constraints"].append({
                    "type": "negation",
                    "details": {
                        "target": "An intermediate molecule lacks a halogen functional group."
                    }
                })
                return

        # Traverse children
        for child in node.get("children", []):
            # New logic for depth calculation
            new_depth = depth
            if node["type"] != "reaction": # If current node is not 'reaction' (e.g., 'mol'), increase depth
                new_depth = depth + 1
            
            dfs_traverse(child, new_depth)

    # Start traversal and set initial value
    halogen_maintained = True
    dfs_traverse(route)

    if halogen_maintained:
        print("Halogen substituent maintained throughout synthesis")
    else:
        print("Halogen substituent NOT maintained throughout synthesis")

    return halogen_maintained, findings_json