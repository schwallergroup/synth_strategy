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
    This function detects late-stage alcohol activation via mesylation.
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

    found_mesylation = False

    def dfs_traverse(node, depth=0):
        nonlocal found_mesylation, findings_json

        if (
            node.get("type") == "reaction" and depth <= 1
        ):
            # The original implementation using raw RDKit was a source of false positives
            # and has been replaced with robust checker functions that verify the transformation.
            if checker.functional_group_lost(node, "alcohol"):
                findings_json["atomic_checks"]["functional_groups"].append("alcohol")
            if checker.functional_group_formed(node, "mesylate"):
                findings_json["atomic_checks"]["functional_groups"].append("mesylate")

            if checker.functional_group_lost(node, "alcohol") and checker.functional_group_formed(node, "mesylate"):
                found_mesylation = True
                findings_json["atomic_checks"]["named_reactions"].append("alcohol_mesylation")
                # Assuming 'depth <= 1' implies 'last_two_stages' for the positional constraint
                findings_json["structural_constraints"].append({
                    "type": "positional",
                    "details": {
                        "target": "alcohol_mesylation",
                        "position": "last_two_stages"
                    }
                })

        for child in node.get("children", []):
            if node.get("type") == "reaction":
                dfs_traverse(child, depth)
            else:
                dfs_traverse(child, depth + 1)

    dfs_traverse(route)
    return found_mesylation, findings_json