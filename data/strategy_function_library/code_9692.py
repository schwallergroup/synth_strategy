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
    This function detects if the synthesis incorporates a trifluoromethyl-containing
    building block early in the synthesis.
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

    trifluoromethyl_early = False
    max_depth = 0

    # First pass to find maximum depth and assign depth to each node
    def find_max_depth(node, current_depth=0):
        nonlocal max_depth
        max_depth = max(max_depth, current_depth)

        # Assign depth to the current node
        node["depth"] = current_depth

        for child in node.get("children", []):
            # New logic for depth calculation:
            # Depth increases only when traversing from a chemical node to a reaction node.
            # Depth remains the same when traversing from a reaction node to a chemical node.
            if node["type"] == "reaction":
                # If current node is a reaction, depth remains the same for its children (chemicals)
                find_max_depth(child, current_depth)
            else:
                # If current node is a chemical, depth increases for its children (reactions)
                find_max_depth(child, current_depth + 1)

    # Create a deep copy of the route to avoid modifying the original
    route_copy = copy.deepcopy(route)
    find_max_depth(route_copy)

    # Check if the final product contains a trifluoromethyl group
    final_product = route_copy["smiles"]
    final_product_has_cf3 = checker.check_fg("Trifluoro group", final_product)

    if final_product_has_cf3:
        findings_json["atomic_checks"]["functional_groups"].append("Trifluoro group")
        findings_json["structural_constraints"].append({
            "type": "positional",
            "details": {
                "target": "Trifluoro group",
                "position": "final_product"
            }
        })
    else:
        return False, findings_json

    # Second pass to check for trifluoromethyl groups in molecules
    def dfs_traverse(node):
        nonlocal trifluoromethyl_early, findings_json

        if node["type"] == "mol":
            # Check if this molecule contains a trifluoromethyl group
            if checker.check_fg("Trifluoro group", node["smiles"]):
                # Get the depth of this node
                depth = node.get("depth", 0)
                is_building_block = node.get("in_stock", False)

                if is_building_block:
                    findings_json["atomic_checks"]["functional_groups"].append("Trifluoro group")
                    # Consider it "early" if it's in the deeper half of the synthesis
                    # Early in synthesis = high depth value
                    if depth > max_depth / 2:
                        trifluoromethyl_early = True
                        findings_json["structural_constraints"].append({
                            "type": "positional",
                            "details": {
                                "target": "Trifluoro group",
                                "on": "building_block",
                                "position": "first_half"
                            }
                        })

        # If we've already found what we're looking for, we can stop early.
        if trifluoromethyl_early:
            return

        for child in node.get("children", []):
            dfs_traverse(child)

    dfs_traverse(route_copy)

    return trifluoromethyl_early, findings_json
