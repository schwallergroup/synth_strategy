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
    This function detects a synthetic strategy involving both ring opening and
    ring formation steps in the same route.
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

    has_ring_opening = False
    has_ring_formation = False

    # Placeholder for checker object, assuming it's defined elsewhere or passed in a real scenario
    class Checker:
        def is_ring_opening(self, reaction):
            # Dummy implementation
            return False

        def is_ring_formation(self, reaction):
            # Dummy implementation
            return False

    checker = Checker()

    def dfs_traverse(node, depth):
        nonlocal has_ring_opening, has_ring_formation, findings_json

        if node["type"] == "reaction":
            reaction = node.get("reaction")
            if not reaction:
                return

            # Check for ring opening or formation using robust checkers
            if checker.is_ring_opening(reaction):
                has_ring_opening = True
                findings_json["atomic_checks"]["named_reactions"].append("ring_destruction")
            elif checker.is_ring_formation(reaction):
                has_ring_formation = True
                findings_json["atomic_checks"]["named_reactions"].append("ring_formation")

        # Traverse children
        for child in node.get("children", []):
            # New logic for depth calculation
            if node["type"] == "reaction":
                # If current node is a reaction, depth remains the same for children
                dfs_traverse(child, depth)
            else:
                # If current node is not a reaction (e.g., chemical), depth increases
                dfs_traverse(child, depth + 1)

    # Start traversal with initial depth 0
    dfs_traverse(route, 0)

    result = has_ring_opening and has_ring_formation

    if result:
        findings_json["structural_constraints"].append({
            "type": "co-occurrence",
            "details": {
                "targets": [
                    "ring_destruction",
                    "ring_formation"
                ]
            }
        })

    return result, findings_json
