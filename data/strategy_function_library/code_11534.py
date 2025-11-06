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
from steerable_retro.utils.check import Check
from steerable_retro.utils import fuzzy_dict, check

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


SUZUKI_REACTION_TYPES = [
    "Suzuki coupling with boronic acids",
    "Suzuki coupling with boronic esters",
    "Suzuki coupling with boronic acids OTf",
    "Suzuki coupling with boronic esters OTf",
    "Suzuki coupling with sulfonic esters",
    "Suzuki",
]

def main(route) -> Tuple[bool, Dict]:
    """
    Detects a late-stage Suzuki coupling reaction. A reaction is flagged if it
    occurs within the first two steps of the synthesis (depth <= 2) and is
    identified by `checker.check_reaction` as one of the types listed in
    `SUZUKI_REACTION_TYPES`.
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

    has_suzuki_coupling = False

    def dfs_traverse(node, current_depth=0):
        nonlocal has_suzuki_coupling, findings_json
        if has_suzuki_coupling:
            return

        if node["type"] == "reaction" and "metadata" in node and "mapped_reaction_smiles" in node["metadata"]:
            # Per the problem description, depth=1 is the final step.
            # The traversal begins at the root (final product) with current_depth=0.
            # The first reaction(s) are at current_depth=1.
            # This checks for reactions in the final two stages of the synthesis.
            if current_depth <= 2:
                rsmi = node["metadata"]["mapped_reaction_smiles"]
                for suzuki_type in SUZUKI_REACTION_TYPES:
                    if checker.check_reaction(suzuki_type, rsmi):
                        has_suzuki_coupling = True
                        findings_json["atomic_checks"]["named_reactions"].append(suzuki_type)
                        # Add the structural constraint if the condition is met
                        findings_json["structural_constraints"].append({
                            "type": "positional",
                            "details": {
                                "target": "Suzuki coupling",
                                "position": "late_stage",
                                "condition": "The reaction must occur within the last two steps of the synthesis."
                            }
                        })
                        return

        # Continue traversing
        for child in node.get("children", []):
            # New logic: depth increases only when going from chemical to reaction
            # Depth remains the same when going from reaction to chemical
            next_depth = current_depth
            if node["type"] != "reaction": # This means current node is 'chemical'
                next_depth = current_depth + 1
            dfs_traverse(child, next_depth)

    dfs_traverse(route)
    return has_suzuki_coupling, findings_json
