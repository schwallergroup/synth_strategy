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


PRESERVED_HETEROCYCLES = ["indole", "furan"]

def main(route) -> Tuple[bool, Dict]:
    """
    This function detects if the synthesis route preserves specific heterocyclic motifs, as defined in the `PRESERVED_HETEROCYCLES` list. A motif is considered preserved if it is present in molecules at two or more different depths of the synthesis tree.
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

    # Track heterocycles through the synthesis using a dictionary
    heterocycle_depths = {name: [] for name in PRESERVED_HETEROCYCLES}

    def dfs_traverse(node, depth=0):
        if node["type"] == "mol" and "smiles" in node:
            mol_smiles = node["smiles"]

            # Check for each specified heterocycle
            for name in PRESERVED_HETEROCYCLES:
                if checker.check_ring(name, mol_smiles):
                    heterocycle_depths[name].append(depth)
                    if name not in findings_json["atomic_checks"]["ring_systems"]:
                        findings_json["atomic_checks"]["ring_systems"].append(name)

        # Traverse children
        for child in node.get("children", []):
            # New logic for depth calculation
            new_depth = depth
            if node["type"] != "reaction": # Only increase depth if current node is not 'reaction'
                new_depth = depth + 1
            dfs_traverse(child, new_depth)

    # Start traversal from root
    dfs_traverse(route)

    # Check if any heterocycle is preserved throughout the synthesis
    # This means it appears at multiple depths (at least 2 different depths)
    is_any_preserved = any(len(set(depths)) >= 2 for depths in heterocycle_depths.values())

    if is_any_preserved:
        findings_json["structural_constraints"].append({
            "type": "count",
            "details": {
                "target": "unique_depths_for_any_target_ring",
                "operator": ">=",
                "value": 2
            }
        })

    return is_any_preserved, findings_json
