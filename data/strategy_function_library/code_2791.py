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
    This function detects early introduction of a trifluoro group
    that is carried through the synthesis.
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

    max_depth = 0
    trifluoroethyl_appearances = {}  # Track depths where trifluoroethyl groups appear
    result = False

    def dfs_traverse(node, depth=0):
        nonlocal max_depth, trifluoroethyl_appearances, findings_json
        max_depth = max(max_depth, depth)

        if node["type"] == "mol" and "smiles" in node:
            # Check for trifluoro group using the checker function
            if checker.check_fg("Trifluoro group", node["smiles"]):
                print(f"Found trifluoro group at depth {depth} in molecule: {node['smiles']}")
                trifluoroethyl_appearances[depth] = node["smiles"]
                if "Trifluoro group" not in findings_json["atomic_checks"]["functional_groups"]:
                    findings_json["atomic_checks"]["functional_groups"].append("Trifluoro group")

        # Continue traversing the tree
        for child in node.get("children", []):
            # New logic for depth calculation
            new_depth = depth
            if node["type"] != "reaction": # If current node is not 'reaction' (e.g., 'chemical' or 'mol')
                new_depth = depth + 1
            dfs_traverse(child, new_depth)

    # Start traversal from the root
    dfs_traverse(route)

    # Check if trifluoroethyl was introduced early and carried through
    if trifluoroethyl_appearances and max_depth > 0:
        # Early introduction means it appears in the deeper half of the synthesis
        # (higher depth values = earlier in synthesis in retrosynthetic analysis)
        early_threshold = max_depth / 2

        # Find the earliest appearance of the trifluoroethyl group
        earliest_appearance = max(trifluoroethyl_appearances.keys())

        # Check if it appears early in the synthesis (deeper than threshold)
        early_introduction = earliest_appearance > early_threshold

        # Check if it persists through to the final product (depth 0)
        persists_to_final = 0 in trifluoroethyl_appearances

        print(f"Earliest trifluoroethyl appearance: depth {earliest_appearance}")
        print(f"Early threshold: {early_threshold}")
        print(f"Early introduction: {early_introduction}")
        print(f"Persists to final product: {persists_to_final}")

        # Update findings_json based on structural constraints
        if early_introduction:
            findings_json["structural_constraints"].append({
                "type": "positional",
                "details": {
                    "target": "Trifluoro group",
                    "position": "early_stage_introduction"
                }
            })
        if persists_to_final:
            findings_json["structural_constraints"].append({
                "type": "positional",
                "details": {
                    "target": "Trifluoro group",
                    "position": "last_stage"
                }
            })

        # Set the main result boolean
        result = early_introduction and persists_to_final

    else:
        print("No trifluoroethyl group found or synthesis has zero depth")

    return result, findings_json
