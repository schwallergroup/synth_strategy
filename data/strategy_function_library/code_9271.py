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


NITROGEN_FGS_OF_INTEREST = ['sulfonamide', 'amide', 'carbamate', 'amine']

def main(route) -> Tuple[bool, Dict]:
    """
    Detects a synthetic strategy involving a sequence of at least two reactions that form or cleave specific nitrogen-containing functional groups. The groups checked are defined in the `NITROGEN_FGS_OF_INTEREST` list, which includes: sulfonamide, amide, carbamate, and amine.
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

    # Track nitrogen modifications
    n_modifications = 0
    n_modification_depths = []
    result = False

    def dfs_traverse(node, depth=0):
        nonlocal n_modifications, n_modification_depths, findings_json

        if node["type"] == "reaction":
            # Assume the reaction object is available on the node as per standard data structure
            reaction = node["reaction"]

            # Corrected logic using checker API and the enumerated list.
            # This replaces the entire buggy, SMARTS-based try-except block.
            n_modified = False
            for fg in NITROGEN_FGS_OF_INTEREST:
                if checker.is_fg_formed(reaction, fg):
                    n_modified = True
                    findings_json["atomic_checks"]["functional_groups"].append(fg)
                    findings_json["atomic_checks"]["named_reactions"].append("functional_group_formation")
                if checker.is_fg_cleaved(reaction, fg):
                    n_modified = True
                    findings_json["atomic_checks"]["functional_groups"].append(fg)
                    findings_json["atomic_checks"]["named_reactions"].append("functional_group_cleavage")

            if n_modified:
                n_modifications += 1
                n_modification_depths.append(depth)

        # Process children
        for child in node.get("children", []):
            # New depth calculation logic
            if node["type"] == "reaction":
                # If current node is a reaction, depth remains the same for children
                dfs_traverse(child, depth)
            else:
                # If current node is not a reaction (e.g., chemical), depth increases
                dfs_traverse(child, depth + 1)

    # Start traversal
    dfs_traverse(route)

    # Check if the strategy is present
    # We need at least 2 nitrogen modifications at different depths
    if n_modifications >= 2:
        findings_json["structural_constraints"].append({
            "type": "count",
            "details": {
                "target": "modification_of_specified_nitrogen_fg",
                "operator": ">=",
                "value": 2
            }
        })
        if len(set(n_modification_depths)) >= 2:
            findings_json["structural_constraints"].append({
                "type": "count",
                "details": {
                    "target": "distinct_stages_with_nitrogen_fg_modification",
                    "operator": ">=",
                    "value": 2
                }
            })
            result = True

    return result, findings_json
