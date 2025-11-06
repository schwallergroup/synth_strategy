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
    This function detects a synthetic strategy involving a morpholine scaffold
    that is maintained throughout the synthesis.
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

    # Track if morpholine is present in the final product and intermediates
    morpholine_in_final = False
    morpholine_intermediates = []

    def dfs_traverse(node, depth=0):
        nonlocal morpholine_in_final, findings_json

        if node["type"] == "mol":
            # Skip empty SMILES
            if not node["smiles"]:
                return

            try:
                # Check if molecule contains morpholine using the checker function
                has_morpholine = checker.check_ring("morpholine", node["smiles"])

                if has_morpholine:
                    findings_json["atomic_checks"]["ring_systems"].append("morpholine")
                    # Check if this is the final product (root node)
                    if depth == 0:
                        morpholine_in_final = True
                        print(f"Found morpholine in final product: {node['smiles']}")
                    # Check if this is an intermediate (not a starting material)
                    elif not node.get("in_stock", False):
                        morpholine_intermediates.append((node["smiles"], depth))
                        print(
                            f"Found morpholine in intermediate at depth {depth}: {node['smiles']}"
                        )
            except Exception as e:
                print(f"Error processing SMILES: {node['smiles']} - {str(e)}")

        # Determine the new depth for recursive calls
        new_depth = depth
        if node["type"] != "reaction":  # If current node is 'mol' (chemical)
            new_depth = depth + 1

        # Traverse children
        for child in node.get("children", []):
            dfs_traverse(child, new_depth)

    # Start traversal
    dfs_traverse(route)

    # Strategy is present if morpholine is in both final product and at least two intermediates
    strategy_present = morpholine_in_final and len(morpholine_intermediates) >= 2

    if morpholine_in_final:
        findings_json["structural_constraints"].append({
            "type": "positional",
            "details": {
                "target": "morpholine",
                "position": "last_stage"
            }
        })
    if len(morpholine_intermediates) >= 2:
        findings_json["structural_constraints"].append({
            "type": "count",
            "details": {
                "target": "morpholine_in_intermediate",
                "operator": ">=",
                "value": 2
            }
        })

    print(f"Morpholine scaffold strategy detected: {strategy_present}")
    print(f"Found {len(morpholine_intermediates)} intermediates with morpholine")

    return strategy_present, findings_json
