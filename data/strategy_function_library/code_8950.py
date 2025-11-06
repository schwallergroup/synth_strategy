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
    Detects if the synthesis route involves early-stage fragment coupling via ether linkage.
    Early stage means at high depth in the tree (further from the final product).
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

    ether_formation_found = False
    ether_formation_depth = -1

    def dfs_traverse(node, depth=0):
        nonlocal ether_formation_found, ether_formation_depth, findings_json

        if node["type"] == "reaction":
            reaction = node.get("reaction")
            # Use a robust checker for any ether formation, which is more accurate.
            if reaction and checker.is_fg_formed(reaction, "ether"):
                ether_formation_found = True
                ether_formation_depth = max(ether_formation_depth, depth)
                findings_json["atomic_checks"]["functional_groups"].append("ether")
                findings_json["atomic_checks"]["named_reactions"].append("ether_formation")

        # Traverse children
        for child in node.get("children", []):
            # New logic for depth calculation
            if node["type"] == "reaction":
                dfs_traverse(child, depth)
            else:
                dfs_traverse(child, depth + 1)

    dfs_traverse(route)

    # Consider it early-stage if it happens at depth 4 or higher
    result = ether_formation_found and ether_formation_depth >= 4

    if result:
        findings_json["structural_constraints"].append({
            "type": "positional",
            "details": {
                "target": "ether_formation",
                "min_depth": 4
            }
        })

    return result, findings_json
