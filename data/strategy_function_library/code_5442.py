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

def main(route) -> Tuple[bool, Dict]:
    """
    Detects if the synthesis route uses an early-stage esterification
    (carboxylic acid to ethyl ester transformation).
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

    found_esterification = False

    # In a real implementation, max_depth would be pre-calculated from the route.
    # We assume a hypothetical function get_max_depth(route) exists.
    max_depth = get_max_depth(route) if route else 0

    def dfs_traverse(node, depth, max_depth):
        nonlocal found_esterification, findings_json

        # The check `depth >= 2` correctly identifies any step that is not the final step (depth=1).
        if node.get("type") == "reaction" and depth >= 2:
            # Correctly check for the specific transformation using the checker API.
            # This avoids the false positives of the original implementation.
            if checker.is_fg_consumed(node, "carboxylic_acid") and checker.is_fg_formed(
                node, "ethyl_ester"
            ):
                found_esterification = True
                findings_json["atomic_checks"]["functional_groups"].append("carboxylic_acid")
                findings_json["atomic_checks"]["functional_groups"].append("ethyl_ester")
                findings_json["atomic_checks"]["named_reactions"].append("esterification")
                # Add the structural constraint if the esterification is found at an early stage
                findings_json["structural_constraints"].append({
                    "type": "positional",
                    "details": {
                        "target": "esterification",
                        "position": "not_last_stage"
                    }
                })

        # Traverse children, propagating context correctly.
        for child in node.get("children", []):
            # New logic for depth calculation:
            # Depth increases only when traversing from a chemical node to a reaction node.
            # Depth remains the same when traversing from a reaction node to a chemical node.
            if node.get("type") == "reaction":
                dfs_traverse(child, depth, max_depth)
            else:  # Assuming 'chemical' or other non-reaction type
                dfs_traverse(child, depth + 1, max_depth)

    # Start traversal from the root node (target molecule) at depth=0.
    if route:
        dfs_traverse(route, 0, max_depth)
        
    return found_esterification, findings_json
