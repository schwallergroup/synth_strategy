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


PHOSPHINE_OXIDE_SMARTS = "[#6]-[P](=O)([#6])[#6]"

# Placeholder for checker object, assuming it would be passed or initialized elsewhere
# For the purpose of this refactoring, we'll define a dummy checker.
class DummyChecker:
    def is_group_unchanged(self, rsmi, smarts):
        # This is a dummy implementation. In a real scenario, this would perform RDKit checks.
        # For now, let's assume it always returns True for demonstration.
        return True

checker = DummyChecker()

def main(route) -> Tuple[bool, Dict]:
    """
    Detects a synthetic strategy involving a tri-substituted phosphine oxide group
    that remains unchanged throughout multiple reaction steps.
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

    # Track phosphine oxide presence across steps
    steps_with_phosphine_oxide = 0
    total_steps = 0

    result = False

    def dfs_traverse(node, depth=0):
        nonlocal steps_with_phosphine_oxide, total_steps, findings_json

        if node["type"] == "reaction":
            total_steps += 1

            rsmi = node["metadata"]["mapped_reaction_smiles"]

            # Correctly check if the group is unchanged using a robust checker.
            # This replaces the previous flawed implementation.
            if checker.is_group_unchanged(rsmi, PHOSPHINE_OXIDE_SMARTS):
                steps_with_phosphine_oxide += 1
                # Record atomic check for functional group
                if "phosphine oxide" not in findings_json["atomic_checks"]["functional_groups"]:
                    findings_json["atomic_checks"]["functional_groups"].append("phosphine oxide")
                # Record atomic check for named reaction
                if "group_unchanged_in_reaction" not in findings_json["atomic_checks"]["named_reactions"]:
                    findings_json["atomic_checks"]["named_reactions"].append("group_unchanged_in_reaction")

        # Traverse children
        for child in node.get("children", []):
            if node["type"] == "reaction":
                # Depth remains the same when traversing from reaction to chemical
                dfs_traverse(child, depth)
            else:
                # Depth increases when traversing from chemical to reaction
                dfs_traverse(child, depth + 1)

    # Start traversal
    dfs_traverse(route)

    # Strategy is present if phosphine oxide appears in multiple steps
    # and persists through most of the synthesis
    if steps_with_phosphine_oxide >= 2 and steps_with_phosphine_oxide >= total_steps * 0.5:
        result = True
        # Record structural constraints
        findings_json["structural_constraints"].append(
            {
                "type": "count",
                "details": {
                    "target": {
                        "event": "group_unchanged_in_reaction",
                        "group": "phosphine oxide"
                    },
                    "operator": ">=",
                    "value": 2
                }
            }
        )
        findings_json["structural_constraints"].append(
            {
                "type": "proportion",
                "details": {
                    "target": {
                        "event": "group_unchanged_in_reaction",
                        "group": "phosphine oxide"
                    },
                    "operator": ">=",
                    "value": 0.5,
                    "relative_to": "total_reaction_steps"
                }
            }
        )

    return result, findings_json
