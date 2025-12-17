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


SUZUKI_REACTION_NAMES = [
    "Suzuki coupling with boronic acids",
    "Suzuki coupling with boronic esters",
    "Suzuki coupling with boronic acids OTf",
    "Suzuki coupling with boronic esters OTf",
    "Suzuki coupling with sulfonic esters",
    "{Suzuki}",
    "Suzuki",
]

def main(route) -> Tuple[bool, Dict]:
    """
    Detects a late-stage Suzuki coupling. This is achieved by checking if any reaction step
    within the final three steps of the synthesis (depth <= 3, where depth=1 is the final step)
    is identified as a Suzuki reaction. The specific reaction names checked are defined in the
    `SUZUKI_REACTION_NAMES` list.
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

    suzuki_coupling_found = False
    suzuki_depth = float("inf")

    def dfs_traverse(node, depth=0):
        nonlocal suzuki_coupling_found, suzuki_depth, findings_json

        if node["type"] == "reaction":
            # The first reaction node is at depth=1
            reaction_depth = depth

            if "mapped_reaction_smiles" in node.get("metadata", {}):
                rsmi = node["metadata"]["mapped_reaction_smiles"]

                # Check for Suzuki coupling using the checker function against a defined list
                is_suzuki = False
                for name in SUZUKI_REACTION_NAMES:
                    if checker.check_reaction(name, rsmi):
                        is_suzuki = True
                        findings_json["atomic_checks"]["named_reactions"].append(name)
                        break

                if is_suzuki:
                    suzuki_coupling_found = True
                    suzuki_depth = min(suzuki_depth, reaction_depth)

        for child in node.get("children", []):
            # New logic: depth increases only when traversing from chemical to reaction
            new_depth = depth
            if node["type"] != "reaction": # This means current node is 'chemical'
                new_depth = depth + 1
            dfs_traverse(child, new_depth)

    # The first call to dfs_traverse is on the root node (final product), with depth=0.
    # The first reaction node(s) will be at depth=1.
    dfs_traverse(route)

    # Consider it late-stage if it occurs at depth 1, 2, or 3
    result = suzuki_coupling_found and suzuki_depth <= 3

    if result:
        # Add the structural constraint if the condition is met
        findings_json["structural_constraints"].append({
            "type": "positional",
            "details": {
                "targets": [
                    "Suzuki coupling with boronic acids",
                    "Suzuki coupling with boronic esters",
                    "Suzuki coupling with boronic acids OTf",
                    "Suzuki coupling with boronic esters OTf",
                    "Suzuki coupling with sulfonic esters",
                    "{Suzuki}",
                    "Suzuki"
                ],
                "position": "within_last_n_stages",
                "value": 3
            }
        })

    return result, findings_json
