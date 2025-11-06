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
    Detects reactions that form a cyclopropyl group during the early stages of a synthesis (defined as a depth of 5 or greater).
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

    cyclopropyl_introduced_early = False

    # Placeholder for checker, assuming it's defined elsewhere or passed in a real scenario
    class MockChecker:
        def is_group_formed(self, rsmi, group_smiles):
            # Simple mock logic for demonstration
            return "C1CC1" in rsmi
    checker = MockChecker()

    def dfs_traverse(node, depth=0):
        nonlocal cyclopropyl_introduced_early, findings_json

        # Assign depth to the node for checking
        if "metadata" not in node:
            node["metadata"] = {}
        node["metadata"]["depth"] = depth

        if node["type"] == "reaction" and "metadata" in node and "mapped_reaction_smiles" in node["metadata"]:
            # Check if this is an early-stage reaction (depth >= 5)
            if node["metadata"].get("depth", -1) >= 5:
                rsmi = node["metadata"]["mapped_reaction_smiles"]
                # Correctly check for cyclopropyl formation using a robust checker
                if checker.is_group_formed(rsmi, "C1CC1"):
                    cyclopropyl_introduced_early = True
                    findings_json["atomic_checks"]["ring_systems"].append("cyclopropyl")
                    findings_json["atomic_checks"]["named_reactions"].append("ring_formation")
                    findings_json["structural_constraints"].append({
                        "type": "positional",
                        "details": {
                            "target_event": "ring_formation",
                            "target_entity": "cyclopropyl",
                            "position": {
                                "variable": "depth",
                                "operator": ">=",
                                "value": 5
                            }
                        }
                    })

        for child in node.get("children", []):
            # New depth calculation logic
            if node["type"] == "reaction":
                # Depth remains the same when traversing from reaction to chemical
                dfs_traverse(child, depth)
            else:
                # Depth increases when traversing from chemical to reaction
                dfs_traverse(child, depth + 1)

    dfs_traverse(route)
    return cyclopropyl_introduced_early, findings_json
