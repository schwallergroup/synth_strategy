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
from synth_strategy.utils.check import Check
from synth_strategy.utils import fuzzy_dict, check

from pathlib import Path
root_data = Path(__file__).parent.parent

fg_args = {
    "file_path": f"{root_data}/patterns/functional_groups.json",
    "value_field": "pattern",
    "key_field": "name",
}
reaction_class_args = {
    "file_path": f"{root_data}/patterns/smirks.json",
    "value_field": "smirks",
    "key_field": "name",
}
ring_smiles_args = {
    "file_path": f"{root_data}/patterns/chemical_rings_smiles.json",
    "value_field": "smiles",
    "key_field": "name",
}
functional_groups = fuzzy_dict.FuzzyDict.from_json(**fg_args)
reaction_classes = fuzzy_dict.FuzzyDict.from_json(**reaction_class_args)
ring_smiles = fuzzy_dict.FuzzyDict.from_json(**ring_smiles_args)

checker = check.Check(
    fg_dict=functional_groups, reaction_dict=reaction_classes, ring_dict=ring_smiles
)


def main(route) -> Tuple[bool, Dict]:
    """
    This function detects if a trifluoromethyl-containing molecule is present
    in the late stage of the synthesis (low depth).
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

    trifluoromethyl_depth = None
    result = False

    def dfs_traverse(node, depth=0):
        nonlocal trifluoromethyl_depth, findings_json

        if node["type"] == "mol":
            # Check if molecule contains trifluoromethyl group
            mol_smiles = node["smiles"]
            if checker.check_fg("Trifluoro group", mol_smiles):
                findings_json["atomic_checks"]["functional_groups"].append("Trifluoro group")
                # Track the lowest depth (latest stage) where trifluoromethyl appears
                if trifluoromethyl_depth is None or depth < trifluoromethyl_depth:
                    trifluoromethyl_depth = depth

        # Determine the depth for the next recursive call
        next_depth = depth
        if node["type"] != "reaction":  # If current node is 'mol' (chemical), depth increases
            next_depth = depth + 1

        # Traverse children
        for child in node.get("children", []):
            dfs_traverse(child, next_depth)

    # Start traversal
    dfs_traverse(route)

    # Return True if trifluoromethyl group is found in late stage (depth <= 2)
    result = trifluoromethyl_depth is not None and trifluoromethyl_depth <= 2

    if result:
        # Add the structural constraint if the condition is met
        findings_json["structural_constraints"].append({
            "type": "positional",
            "details": {
                "target": "Trifluoro group",
                "position": "late_stage",
                "condition": "The latest introduction of the group must be at a depth <= 2."
            }
        })

    return result, findings_json
