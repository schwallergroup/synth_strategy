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

def main(route) -> Tuple[bool, Dict]:
    """
    This function detects a late-stage Suzuki coupling strategy, defined as a Suzuki reaction
    occurring in the latter half of the synthetic sequence (i.e., at a depth less than or
    equal to half the maximum depth).
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

    suzuki_found = False
    suzuki_depth = -1
    max_depth = -1

    def dfs_traverse(node, depth=0):
        nonlocal suzuki_found, suzuki_depth, max_depth, findings_json

        max_depth = max(max_depth, depth)

        if node["type"] == "reaction":
            # Use the robust checker for named reactions, which is not prone to FPs.
            if checker.is_named_reaction(node, "suzuki"):
                suzuki_found = True
                suzuki_depth = depth
                findings_json["atomic_checks"]["named_reactions"].append("suzuki")
                print(f"Found Suzuki coupling at depth {depth}")

        for child in node.get("children", []):
            # New logic: depth increases only when traversing from chemical to reaction.
            # Depth remains the same when traversing from reaction to chemical.
            if node["type"] == "reaction":
                dfs_traverse(child, depth)
            else:  # Assuming 'chemical' or other types that should increase depth
                dfs_traverse(child, depth + 1)

    dfs_traverse(route)

    # Consider it late-stage if it occurs in the first half of the synthesis (lower depth)
    is_late_stage = suzuki_found and suzuki_depth != -1 and suzuki_depth <= max_depth / 2

    if is_late_stage:
        findings_json["structural_constraints"].append({
            "type": "positional",
            "details": {
                "target": "suzuki",
                "position": "latter_half"
            }
        })
        print(
            f"Detected late-stage Suzuki coupling strategy at depth {suzuki_depth} (max depth: {max_depth})"
        )

    return is_late_stage, findings_json
