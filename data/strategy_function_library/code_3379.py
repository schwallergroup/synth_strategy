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

def main(route) -> Tuple[bool, Dict]:
    """
    Detects syntheses where a morpholine group is present in early-stage intermediates
    but is removed before the final stages, indicating its potential use as a transient group.
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

    # Initialize tracking variables
    has_morpholine = False
    morpholine_depth = float("inf")  # Initialize to infinity
    max_depth = 0

    def dfs_traverse(node, depth=0):
        nonlocal has_morpholine, morpholine_depth, max_depth, findings_json

        # Track maximum depth
        max_depth = max(max_depth, depth)

        if node["type"] == "mol":
            # Check for morpholine in molecule using the checker function
            if checker.check_ring("morpholine", node["smiles"]):
                has_morpholine = True
                findings_json["atomic_checks"]["ring_systems"].append("morpholine")
                # Update morpholine_depth to the minimum depth where it's found
                morpholine_depth = min(morpholine_depth, depth)

        # Continue traversing
        for child in node.get("children", []):
            # New logic for depth calculation
            if node["type"] == "reaction":
                dfs_traverse(child, depth)
            else:  # node['type'] is 'mol' or 'chemical'
                dfs_traverse(child, depth + 1)

    # Start traversal
    dfs_traverse(route)

    # If morpholine was never found, set depth to 0
    if morpholine_depth == float("inf"):
        morpholine_depth = 0

    # Strategy is present if the latest appearance of morpholine (min_depth)
    # is in the first half of the synthesis (depth > max_depth / 2).
    # This implies it was removed before the final stages.
    early_incorporation = has_morpholine and morpholine_depth > (max_depth / 2)

    if has_morpholine:
        findings_json["structural_constraints"].append({
            "type": "co-occurrence",
            "details": {
                "targets": [
                    "morpholine"
                ]
            }
        })

    if early_incorporation:
        findings_json["structural_constraints"].append({
            "type": "positional",
            "details": {
                "target": "morpholine",
                "position": "early_half"
            }
        })

    return early_incorporation, findings_json
