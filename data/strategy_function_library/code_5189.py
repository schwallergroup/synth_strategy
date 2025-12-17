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

HETEROCYCLES_TO_PRESERVE = ['benzimidazole', 'quinoline']

def main(route) -> Tuple[bool, Dict]:
    """
    This function detects if the synthesis route preserves specific heterocyclic structures
    (defined in HETEROCYCLES_TO_PRESERVE) throughout the synthesis.
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
    result = True

    # Track molecules at each depth
    molecules_by_depth = {}

    # Add depth information during traversal
    def add_depth_info(node, current_depth=0):
        node["depth"] = current_depth
        for child in node.get("children", []):
            # New logic for depth calculation
            if node["type"] == "reaction":
                # Depth remains the same when traversing from reaction to chemical
                add_depth_info(child, current_depth)
            else:
                # Depth increases when traversing from chemical to reaction
                add_depth_info(child, current_depth + 1)

    # Add depth information to the route
    add_depth_info(route)

    def dfs_traverse(node):
        nonlocal findings_json
        if node["type"] == "mol":
            depth = node.get("depth", 0)

            if depth not in molecules_by_depth:
                molecules_by_depth[depth] = []
            molecules_by_depth[depth].append(node["smiles"])

        # Traverse children
        for child in node.get("children", []):
            dfs_traverse(child)

    # Start traversal from the root
    dfs_traverse(route)

    # Get depths in order (0 is the target molecule)
    depths = sorted(molecules_by_depth.keys())
    if not depths:
        return False, findings_json

    # Check which heterocycles are in the target molecule
    target_depth = depths[0]
    target_has_heterocycles = {h: False for h in HETEROCYCLES_TO_PRESERVE}

    for smiles in molecules_by_depth[target_depth]:
        for h in HETEROCYCLES_TO_PRESERVE:
            if not target_has_heterocycles[h] and checker.check_ring(h, smiles):
                target_has_heterocycles[h] = True
                if h not in findings_json["atomic_checks"]["ring_systems"]:
                    findings_json["atomic_checks"]["ring_systems"].append(h)

    # If target doesn't have any of the specified heterocycles, the strategy is not applicable
    if not any(target_has_heterocycles.values()):
        result = False
    else:
        # Add positional constraint if target has at least one of the heterocycles
        findings_json["structural_constraints"].append({
            "type": "positional",
            "details": {
                "target": [
                    "benzimidazole",
                    "quinoline"
                ],
                "position": "last_stage",
                "condition": "at_least_one"
            }
        })

    if not result:
        return result, findings_json

    # Check if heterocycles present in the target are preserved at all depths
    for depth in depths:
        depth_has_heterocycles = {h: False for h in HETEROCYCLES_TO_PRESERVE}
        for smiles in molecules_by_depth[depth]:
            for h in HETEROCYCLES_TO_PRESERVE:
                if not depth_has_heterocycles[h] and checker.check_ring(h, smiles):
                    depth_has_heterocycles[h] = True
                    if h not in findings_json["atomic_checks"]["ring_systems"]:
                        findings_json["atomic_checks"]["ring_systems"].append(h)

        # Check for preservation of heterocycles that are in the target
        for h in HETEROCYCLES_TO_PRESERVE:
            if target_has_heterocycles[h] and not depth_has_heterocycles[h]:
                result = False
                # If a heterocycle is not preserved, it implies a ring destruction event
                if "ring_destruction" not in findings_json["atomic_checks"]["named_reactions"]:
                    findings_json["atomic_checks"]["named_reactions"].append("ring_destruction")
                findings_json["structural_constraints"].append({
                    "type": "negation",
                    "details": {
                        "event": "ring_destruction",
                        "scope": [
                            "benzimidazole",
                            "quinoline"
                        ]
                    }
                })
                return result, findings_json

    return result, findings_json
