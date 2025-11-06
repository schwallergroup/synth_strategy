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

def detect_late_stage_boc_deprotection(reaction, depth, max_depth) -> Tuple[bool, Dict]:
    """
    Checks for the deprotection of a Boc group in the final (depth=1) or penultimate (depth=2) reaction step.
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
    result = False

    # A reaction's depth must be 1 (final) or 2 (penultimate) to be considered late-stage.
    if depth not in [1, 2]:
        return result, findings_json

    # Use a robust checker for Boc deprotection. This is superior to a specific SMARTS pattern
    # as it correctly identifies the removal of any Boc protecting group (e.g., N-Boc, O-Boc).
    is_boc_deprotection = checker.check_functional_group_disappearance(reaction, "Boc")

    if is_boc_deprotection:
        result = True
        findings_json["atomic_checks"]["named_reactions"].append("Boc_deprotection")
        findings_json["atomic_checks"]["functional_groups"].append("Boc")

        # Record structural constraint based on depth
        if depth == 1:
            findings_json["structural_constraints"].append({
                "type": "positional",
                "details": {
                    "target": "Boc_deprotection",
                    "position": [
                        "final_step"
                    ]
                }
            })
        elif depth == 2:
            findings_json["structural_constraints"].append({
                "type": "positional",
                "details": {
                    "target": "Boc_deprotection",
                    "position": [
                        "penultimate_step"
                    ]
                }
            })

    return result, findings_json