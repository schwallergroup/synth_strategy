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
    Detects if the final reaction step involves the formation of an amide functional group.
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

    # The strategy only concerns the final reaction (the root of the synthesis tree).
    # No traversal is necessary; a direct check is sufficient and correct.
    if route and route.get('type') == 'reaction':
        # The final reaction is the root node, equivalent to depth=1.
        if checker.check_functional_group_formation(route, 'amide'):
            result = True
            findings_json["atomic_checks"]["functional_groups"].append("amide")
            findings_json["structural_constraints"].append({
                "type": "positional",
                "details": {
                    "target": "amide_formation",
                    "position": "last_stage"
                }
            })

    return result, findings_json