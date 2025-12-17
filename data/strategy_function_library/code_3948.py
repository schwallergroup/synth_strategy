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
    Detects whether an indole scaffold is preserved throughout the synthesis.

    Returns True if:
    1. The target molecule contains an indole scaffold
    2. The majority of late-stage intermediates preserve the indole scaffold
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

    molecules_with_indole = 0
    total_molecules = 0
    target_has_indole = False
    late_stage_with_indole = 0
    late_stage_total = 0

    # Define what counts as "late stage" (depth < 3)
    LATE_STAGE_DEPTH = 3

    def dfs_traverse(node, depth=0):
        nonlocal molecules_with_indole, total_molecules, target_has_indole
        nonlocal late_stage_with_indole, late_stage_total, findings_json

        if node["type"] == "mol" and node["smiles"]:
            has_indole = checker.check_ring("indole", node["smiles"])
            if has_indole:
                findings_json["atomic_checks"]["ring_systems"].append("indole")

            # Check if this is the target molecule (root node)
            if depth == 0:
                target_has_indole = has_indole
                print(f"Target molecule has indole: {target_has_indole}")

            # Track late-stage intermediates separately
            if depth < LATE_STAGE_DEPTH:
                late_stage_total += 1
                if has_indole:
                    late_stage_with_indole += 1

        for child in node.get("children", []):
            # New logic for depth calculation
            new_depth = depth
            if node["type"] != "reaction": # If current node is chemical, depth increases
                new_depth = depth + 1
            # If current node is reaction, depth remains the same
            dfs_traverse(child, new_depth)

    # Start traversal
    dfs_traverse(route)

    # Calculate preservation ratios
    overall_ratio = molecules_with_indole / total_molecules if total_molecules > 0 else 0
    late_stage_ratio = late_stage_with_indole / late_stage_total if late_stage_total > 0 else 0

    print(f"Molecules with indole: {molecules_with_indole}/{total_molecules} ({overall_ratio:.2f})")
    print(
        f"Late-stage molecules with indole: {late_stage_with_indole}/{late_stage_total} ({late_stage_ratio:.2f})"
    )

    # Scaffold is preserved if:
    # 1. Target molecule has indole AND
    # 2. Either majority of all molecules have indole OR majority of late-stage intermediates have indole
    has_preserved_scaffold = target_has_indole and (overall_ratio >= 0.5 or late_stage_ratio >= 0.5)

    if target_has_indole:
        findings_json["structural_constraints"].append({"type": "positional", "details": {"target": "indole", "position": "last_stage"}})
    
    if late_stage_ratio >= 0.5:
        findings_json["structural_constraints"].append({"type": "count", "details": {"target": "ratio_of_indole_in_late_stage_molecules", "operator": ">=", "value": 0.5}})

    print(f"Indole scaffold preserved: {has_preserved_scaffold}")

    return has_preserved_scaffold, findings_json
