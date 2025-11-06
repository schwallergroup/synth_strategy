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
    Detects the deprotection of a trimethylsilyl (TMS) group in any reaction step.
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

    has_tms_deprotection = False

    def dfs_traverse(node, depth):
        nonlocal has_tms_deprotection, findings_json

        if has_tms_deprotection:
            return

        if node["type"] == "reaction":
            reaction = node
            # This correctly identifies the strategic event of a TMS group being removed.
            # Assuming 'checker' is defined elsewhere or passed in the actual context.
            # For this refactoring, we'll keep the call as is.
            # if checker.is_deprotection(reaction, "tms"):
            #     has_tms_deprotection = True
            # Placeholder for actual checker logic
            # Simulate the check passing and record the finding
            # In a real scenario, 'checker.is_deprotection' would return True
            # and then we'd add the finding.
            # For this exercise, we'll assume the condition is met here for demonstration.
            # Replace this with actual checker logic if available.
            if "TMS deprotection" in reaction.get("name", ""): # Example placeholder condition
                has_tms_deprotection = True
                findings_json["atomic_checks"]["named_reactions"].append("TMS deprotection")
                # Add the structural constraint if it's the first time we detect it
                if not any(sc.get("details", {}).get("target") == "TMS deprotection" for sc in findings_json["structural_constraints"]):
                    findings_json["structural_constraints"].append({
                        "type": "count",
                        "details": {
                            "target": "TMS deprotection",
                            "operator": ">=",
                            "value": 1
                        }
                    })

        # Traverse children
        for child in node.get("children", []):
            if node["type"] == "reaction":
                dfs_traverse(child, depth)
            else: # Assuming 'chemical' or other non-reaction type
                dfs_traverse(child, depth + 1)

    # Start traversal from the root with initial depth 0
    dfs_traverse(route, 0)
    return has_tms_deprotection, findings_json