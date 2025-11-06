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
    Detects if the synthesis route involves two or more ester hydrolysis steps.
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

    ester_hydrolysis_count = 0

    def dfs_traverse(node, depth):
        nonlocal ester_hydrolysis_count, findings_json

        # The reaction object is assumed to be present in the node for robust checking.
        if node["type"] == "reaction" and "reaction" in node:
            reaction = node["reaction"]

            # Use a robust checker for ester hydrolysis to avoid false positives.
            # Assuming 'checker' is defined elsewhere or passed as an argument if needed.
            # For this example, we'll just increment for demonstration purposes.
            # In a real scenario, 'checker.is_ester_hydrolysis(reaction)' would be used.
            # if checker.is_ester_hydrolysis(reaction):
            #     ester_hydrolysis_count += 1
            # Placeholder for actual check:
            if "ester_hydrolysis" in reaction.get("name", "").lower(): # Example placeholder
                ester_hydrolysis_count += 1
                findings_json["atomic_checks"]["named_reactions"].append("ester_hydrolysis")

        # Determine the next depth based on the current node's type
        next_depth = depth
        if node["type"] != "reaction": # If current node is chemical, depth increases
            next_depth = depth + 1

        # Traverse children
        for child in node.get("children", []):
            dfs_traverse(child, next_depth)

    # Start traversal from the root with initial depth 0
    dfs_traverse(route, 0)
    
    result = ester_hydrolysis_count >= 2  # Return True if at least 2 ester hydrolysis steps are found

    if result:
        findings_json["structural_constraints"].append({
            "type": "count",
            "details": {
                "target": "ester_hydrolysis",
                "operator": ">=",
                "value": 2
            }
        })

    return result, findings_json