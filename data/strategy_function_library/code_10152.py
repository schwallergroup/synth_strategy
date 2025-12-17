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
    This function detects if a pyridazine core is maintained throughout the synthesis.
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

    # Track molecules at each depth
    molecules_by_depth = {}

    def dfs_traverse(node, depth=0):
        if node["type"] == "mol":
            smiles = node.get("smiles", "")
            if smiles:
                # Store molecule at this depth
                if depth not in molecules_by_depth:
                    molecules_by_depth[depth] = []
                molecules_by_depth[depth].append(smiles)

        # Determine the next depth based on the current node's type
        next_depth = depth
        if node["type"] != "reaction":  # If current node is 'mol' (chemical), depth increases
            next_depth = depth + 1

        # Traverse children (going deeper in retrosynthesis)
        for child in node.get("children", []):
            dfs_traverse(child, next_depth)

    # Start traversal
    dfs_traverse(route)

    # Check if pyridazine is present at all depths
    has_pyridazine_throughout = True
    all_depths_checked = []
    for depth in sorted(molecules_by_depth.keys()):
        depth_has_pyridazine = False
        for smiles in molecules_by_depth[depth]:
            if checker.check_ring("pyridazine", smiles):
                depth_has_pyridazine = True
                if "pyridazine" not in findings_json["atomic_checks"]["ring_systems"]:
                    findings_json["atomic_checks"]["ring_systems"].append("pyridazine")
                break

        if not depth_has_pyridazine:
            has_pyridazine_throughout = False
            break
        else:
            all_depths_checked.append(depth)

    if has_pyridazine_throughout and all_depths_checked:
        # This structural constraint is met if pyridazine is found at all relevant depths
        findings_json["structural_constraints"].append({
            "type": "positional",
            "details": {
                "target": "pyridazine",
                "position": "all_stages"
            }
        })

    return has_pyridazine_throughout, findings_json
