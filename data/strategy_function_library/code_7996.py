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


SUZUKI_REACTION_VARIANTS = [
    "Suzuki coupling with boronic acids",
    "Suzuki coupling with boronic acids OTf",
    "Suzuki coupling with boronic esters",
    "Suzuki coupling with boronic esters OTf",
    "Suzuki coupling with sulfonic esters",
]

def main(route) -> Tuple[bool, Dict]:
    """
    Checks for a late-stage Suzuki coupling. It identifies if any reaction in the latter half of the synthesis (depth <= max_depth / 2) matches a specific, predefined list of named Suzuki reaction variants.
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

    suzuki_found = False
    max_depth = 0

    # First pass to find max_depth
    def find_max_depth(node, depth=0):
        nonlocal max_depth
        max_depth = max(max_depth, depth)
        for child in node.get("children", []):
            # The depth passed to the recursive call should be `depth` if the current node's type is 'reaction',
            # and `depth + 1` if the current node's type is not 'reaction' (e.g., 'chemical').
            new_depth = depth if node["type"] == "reaction" else depth + 1
            find_max_depth(child, new_depth)

    find_max_depth(route)

    # Second pass to find Suzuki couplings
    def dfs_traverse(node, depth, max_depth):
        nonlocal suzuki_found, findings_json

        if node["type"] == "reaction" and "rsmi" in node.get("metadata", {}):
            rsmi = node["metadata"]["mapped_reaction_smiles"]

            # Check if this is one of the specified Suzuki coupling variants
            is_suzuki_variant_found = False
            for variant in SUZUKI_REACTION_VARIANTS:
                if checker.check_reaction(variant, rsmi):
                    findings_json["atomic_checks"]["named_reactions"].append(variant)
                    is_suzuki_variant_found = True

            if is_suzuki_variant_found:
                # Check if this is in the second half of the synthesis (late stage)
                if depth <= max_depth / 2:
                    suzuki_found = True
                    # Add the structural constraint if a late-stage Suzuki is found
                    findings_json["structural_constraints"].append({
                        "type": "positional",
                        "details": {
                            "target": [
                                "Suzuki coupling with boronic acids",
                                "Suzuki coupling with boronic acids OTf",
                                "Suzuki coupling with boronic esters",
                                "Suzuki coupling with boronic esters OTf",
                                "Suzuki coupling with sulfonic esters"
                            ],
                            "position": "latter_half"
                        }
                    })

        for child in node.get("children", []):
            # The depth passed to the recursive call should be `depth` if the current node's type is 'reaction',
            # and `depth + 1` if the current node's type is not 'reaction' (e.g., 'chemical').
            new_depth = depth if node["type"] == "reaction" else depth + 1
            dfs_traverse(child, new_depth, max_depth)

    dfs_traverse(route, 0, max_depth)
    return suzuki_found, findings_json
