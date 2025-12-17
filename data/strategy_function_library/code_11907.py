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


PRESERVED_HETEROCYCLES_OF_INTEREST = [
    "pyrazole", "pyridine", "benzofuran", "furan", "thiophene",
    "imidazole", "oxazole", "thiazole", "indole", "benzimidazole",
]

def main(route) -> Tuple[bool, Dict]:
    """
    Detects if specific heterocyclic systems from the PRESERVED_HETEROCYCLES_OF_INTEREST list are preserved throughout the synthesis.
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

    # Identify heterocycles in the target molecule (root of the tree)
    target_heterocycles = set()
    if route["type"] == "mol" and "smiles" in route:
        target_smiles = route["smiles"]
        for heterocycle in PRESERVED_HETEROCYCLES_OF_INTEREST:
            if checker.check_ring(heterocycle, target_smiles):
                target_heterocycles.add(heterocycle)
                findings_json["atomic_checks"]["ring_systems"].append(heterocycle)

    if not target_heterocycles:
        return False, findings_json

    # Dictionary to track depths where each heterocycle is found
    heterocycle_depths = {hc: set() for hc in target_heterocycles}

    # Inner function to traverse the synthesis tree and record molecule locations
    def check_preservation(node, depth=0):
        if node["type"] == "mol" and "smiles" in node:
            smiles = node["smiles"]
            for heterocycle in target_heterocycles:
                if checker.check_ring(heterocycle, smiles):
                    heterocycle_depths[heterocycle].add(depth)
                    # Only add to findings_json if it's not already there to avoid duplicates
                    if heterocycle not in findings_json["atomic_checks"]["ring_systems"]:
                        findings_json["atomic_checks"]["ring_systems"].append(heterocycle)

        for child in node.get("children", []):
            # New logic for depth calculation
            new_depth = depth
            if node["type"] != "reaction": # If current node is not 'reaction' (e.g., 'mol'), increase depth
                new_depth = depth + 1
            check_preservation(child, new_depth)

    # Start traversal from the root
    check_preservation(route)

    # The strategy is positive if any target heterocycle appears at more than one depth,
    # indicating it was carried over from a precursor.
    is_preserved = any(len(depths) > 1 for depths in heterocycle_depths.values())

    if is_preserved:
        # Add the structural constraint if the condition is met
        findings_json["structural_constraints"].append({
            "type": "count",
            "details": {
                "target": "depth_count_per_heterocycle",
                "operator": ">",
                "value": 1,
                "scope": "any_target_heterocycle",
                "condition": "The strategy is valid if any heterocycle found in the target molecule is also present in at least one precursor, meaning it appears at more than one depth level in the synthesis tree."
            }
        })

    return is_preserved, findings_json