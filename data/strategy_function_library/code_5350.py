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
    Detects a synthesis strategy centered around a pyridine core that is maintained
    throughout the synthesis.
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

    # Track intermediates with pyridine and total intermediates
    intermediates_with_pyridine = 0
    total_intermediates = 0
    target_has_pyridine = False

    # Track main synthetic path molecules
    main_path_molecules = set()

    def identify_main_path(node, is_main_path=True):
        """Identify molecules on the main synthetic path (containing the core structure)"""
        nonlocal main_path_molecules

        if node["type"] == "mol":
            if is_main_path and not node.get("in_stock", False):
                main_path_molecules.add(node["smiles"])

            # Continue traversing children
            for child in node.get("children", []):
                # If this is a reaction node, all children are on the main path
                if child["type"] == "reaction":
                    identify_main_path(child, is_main_path)

        elif node["type"] == "reaction":
            # For each child (reactant), check if it contributes to the pyridine core
            for child in node.get("children", []):
                if child["type"] == "mol":
                    # If the child has pyridine and the product has pyridine, it's on the main path
                    child_has_pyridine = checker.check_ring("pyridine", child["smiles"])
                    if child_has_pyridine and is_main_path:
                        identify_main_path(child, True)
                    else:
                        # This is likely a reagent, not on the main pyridine path
                        identify_main_path(child, False)
                else:
                    identify_main_path(child, is_main_path)

    def dfs_traverse(node, depth=0):
        nonlocal intermediates_with_pyridine, total_intermediates, target_has_pyridine, findings_json

        if node["type"] == "mol":
            # Target molecule check (depth 0)
            if depth == 0:
                target_has_pyridine = checker.check_ring("pyridine", node["smiles"])
                if target_has_pyridine:
                    findings_json["atomic_checks"]["ring_systems"].append("pyridine")
                    findings_json["structural_constraints"].append({"type": "positional", "details": {"target": "pyridine", "position": "last_stage"}})

            # Only count non-starting materials as intermediates if they're on the main path
            if (
                not node.get("in_stock", False)
                and depth > 0
                and node["smiles"] in main_path_molecules
            ):
                total_intermediates += 1
                # Check if this intermediate has pyridine
                if checker.check_ring("pyridine", node["smiles"]):
                    intermediates_with_pyridine += 1

        # Continue traversing children
        for child in node.get("children", []):
            # New logic for depth calculation
            if node["type"] == "reaction":
                dfs_traverse(child, depth)
            else:  # node['type'] == 'mol'
                dfs_traverse(child, depth + 1)

    # First identify the main synthetic path
    identify_main_path(route)

    # Then traverse and count intermediates
    dfs_traverse(route)

    # Strategy is present if:
    # 1. Target molecule has pyridine
    # 2. All intermediates on the main path maintain pyridine
    # 3. There are intermediates to check
    strategy_present = (
        target_has_pyridine
        and (intermediates_with_pyridine == total_intermediates)
        and total_intermediates > 0
    )

    # Record structural constraints based on final flags
    if intermediates_with_pyridine == total_intermediates:
        findings_json["structural_constraints"].append({"type": "count", "details": {"target": "intermediate_without_pyridine_on_main_path", "operator": "==", "value": 0}})
    
    if total_intermediates > 0:
        findings_json["structural_constraints"].append({"type": "count", "details": {"target": "intermediate", "operator": ">", "value": 0}})

    return strategy_present, findings_json
