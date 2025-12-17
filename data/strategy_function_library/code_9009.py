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
    Detects if the synthetic route employs a late-stage Suzuki coupling. A reaction is considered 'late-stage' if it occurs in the first half of the synthetic sequence (i.e., `depth <= max_depth / 2`).
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

    # Track if we found a Suzuki coupling
    found_suzuki = False
    # Track the depth of the Suzuki coupling
    suzuki_depth = None
    # Track the total depth of the synthesis
    max_depth = 0

    def dfs_traverse(node, depth=0):
        nonlocal found_suzuki, suzuki_depth, max_depth, findings_json

        # Update max depth
        max_depth = max(max_depth, depth)

        # Check if this is a reaction node
        if node.get("type") == "reaction":
            # Check if this is a Suzuki coupling
            if checker.is_named_reaction(node, 'suzuki_coupling'):
                found_suzuki = True
                suzuki_depth = depth
                findings_json["atomic_checks"]["named_reactions"].append("suzuki_coupling")
                findings_json["structural_constraints"].append({"type": "co-occurrence", "details": {"targets": ["suzuki_coupling"]}})
                print(f"Found Suzuki coupling at depth {depth}")

        # Traverse children
        for child in node.get("children", []):
            # New logic: depth increases only when going from chemical to reaction
            # Depth remains the same when going from reaction to chemical
            if node.get("type") == "reaction":
                dfs_traverse(child, depth)
            else: # Assuming 'chemical' or other types that should increase depth
                dfs_traverse(child, depth + 1)

    # Start traversal
    dfs_traverse(route)

    result = False
    # Determine if this is a late-stage Suzuki coupling
    # Late stage means it occurs in the first half of the synthesis (lower depth values)
    if found_suzuki and suzuki_depth is not None:
        # If Suzuki occurs in the first half of the synthesis (considering depth)
        is_late_stage = suzuki_depth <= max_depth / 2
        print(
            f"Suzuki depth: {suzuki_depth}, Max depth: {max_depth}, Is late stage: {is_late_stage}"
        )
        if is_late_stage:
            findings_json["structural_constraints"].append({"type": "positional", "details": {"target": "suzuki_coupling", "position": "first_half_of_route"}})
            result = True

    return result, findings_json