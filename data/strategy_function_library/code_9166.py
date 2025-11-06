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


def dfs_traverse(reaction, depth, max_depth) -> Tuple[bool, Dict]:
    """
    Checks for sulfone formation in the final step (depth=1) of a synthesis.
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

    # A reaction is considered "late-stage" if it's the final step (depth=1).
    if depth == 1:
        # Use the checker API to robustly identify the formation of a sulfone group.
        # This correctly checks for presence in products and absence in reactants,
        # and handles all types of sulfones (C-SO2-C, C-SO2-N, etc.),
        # avoiding the false negatives of the original SMARTS pattern.
        if checker.is_fg_formed(reaction, "sulfone"):
            result = True
            findings_json["atomic_checks"]["functional_groups"].append("sulfone")
            findings_json["atomic_checks"]["named_reactions"].append("functional_group_formation")
            findings_json["structural_constraints"].append({
                "type": "positional",
                "details": {
                    "target": "functional_group_formation",
                    "position": "last_stage"
                }
            })
    return result, findings_json