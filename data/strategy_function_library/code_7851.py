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
    """This function detects if the synthetic route involves a late-stage nitro reduction. A reaction is identified as a nitro reduction using a specific reaction template checker. The strategy is flagged if such a reaction occurs in the final steps of the synthesis (depth < 3, where depth=1 is the final step)."""
    findings_template = {
        "atomic_checks": {
            "named_reactions": [],
            "ring_systems": [],
            "functional_groups": []
        },
        "structural_constraints": []
    }
    findings_json = copy.deepcopy(findings_template)

    nitro_reduction_found = False
    nitro_reduction_depth = float("inf")

    def dfs_traverse(node, depth=0):
        nonlocal nitro_reduction_found, nitro_reduction_depth, findings_json

        if node["type"] == "reaction":
            try:
                # Extract reaction smiles
                rsmi = node["metadata"]["mapped_reaction_smiles"]

                # Check if this is a nitro reduction reaction
                if checker.check_reaction("Reduction of nitro groups to amines", rsmi):
                    print(f"Nitro reduction reaction detected at depth {depth}")
                    nitro_reduction_found = True
                    nitro_reduction_depth = min(nitro_reduction_depth, depth)
                    findings_json["atomic_checks"]["named_reactions"].append("Reduction of nitro groups to amines")
            except Exception as e:
                print(f"Error processing reaction node: {e}")

        # Determine the next depth based on the current node's type
        next_depth = depth
        if node["type"] != "reaction":  # Only increase depth if not a reaction node
            next_depth = depth + 1

        # Traverse children
        for child in node.get("children", []):
            dfs_traverse(child, next_depth)

    # Start traversal
    dfs_traverse(route)

    # Consider it late stage if it's in the first few steps of the synthesis (depth < 3)
    is_late_stage = nitro_reduction_found and nitro_reduction_depth < 3
    print(f"Nitro reduction found: {nitro_reduction_found}, at depth: {nitro_reduction_depth}")
    print(f"Is late stage: {is_late_stage}")

    if is_late_stage:
        # Add the structural constraint if the condition is met
        findings_json["structural_constraints"].append({
            "type": "positional",
            "details": {
                "target": "Reduction of nitro groups to amines",
                "position": {
                    "operator": "<",
                    "value": 3,
                    "unit": "depth_from_product"
                }
            }
        })

    return is_late_stage, findings_json
