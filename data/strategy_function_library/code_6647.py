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


def main(route, checker, max_depth) -> Tuple[bool, Dict]:
    """
    This function detects a synthetic strategy where a diaryl ether motif
    is preserved throughout the synthesis.
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

    diaryl_ether_at_all_depths = True
    depths_checked = set()
    result = False

    def dfs_traverse(node, max_depth, depth=0):
        nonlocal diaryl_ether_at_all_depths, findings_json

        if node["type"] == "reaction":
            # Use robust checker functions instead of manual parsing and specific SMARTS
            product_mol = checker.get_product_mol(node)

            if product_mol:
                has_diaryl_ether = checker.has_functional_group(product_mol, "diaryl ether")
                depths_checked.add(depth)

                if not has_diaryl_ether:
                    diaryl_ether_at_all_depths = False
                    # If diaryl ether is NOT present, this means the 'negation' constraint is met for this step
                    # However, the overall strategy requires it to be present at ALL depths, so this failure
                    # means the 'negation' constraint (absence) is met, which is bad for the strategy.
                    # We record the absence as a finding that contributes to the overall failure.
                    findings_json["structural_constraints"].append({
                        "type": "negation",
                        "details": {
                            "target": "absence_of_diaryl_ether_in_product",
                            "scope": "any_step"
                        }
                    })
                else:
                    # If diaryl ether IS present, record it as an atomic check
                    findings_json["atomic_checks"]["functional_groups"].append("diaryl ether")

        # Traverse children
        for child in node.get("children", []):
            # New logic for depth calculation
            if node["type"] == "reaction":
                # Depth remains the same when traversing from a reaction node to a chemical node
                dfs_traverse(child, max_depth, depth)
            else:
                # Depth increases when traversing from a chemical node to a reaction node
                dfs_traverse(child, max_depth, depth + 1)

    # Start traversal from root. The call is updated to pass max_depth,
    # which is assumed to be available in the calling scope.
    dfs_traverse(route, max_depth)

    # Ensure we checked at least 2 depths
    if len(depths_checked) >= 2 and diaryl_ether_at_all_depths:
        result = True
        # If the strategy is successful, it means the 'count' constraint is met
        findings_json["structural_constraints"].append({
            "type": "count",
            "details": {
                "target": "reaction_steps",
                "operator": ">=",
                "value": 2
            }
        })

    return result, findings_json