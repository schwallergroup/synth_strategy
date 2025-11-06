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


SULFONAMIDE_FORMATION_REACTIONS = [
    "Sulfonamide synthesis (Schotten-Baumann) primary amine",
    "Sulfonamide synthesis (Schotten-Baumann) secondary amine",
]

def main(route) -> Tuple[bool, Dict]:
    """
    Detects if a sulfonamide is formed in the late stages of a synthesis via a specific set of named reactions, including primary and secondary amine-based Schotten-Baumann syntheses.
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

    sulfonamide_found = False
    late_stage = False
    max_depth_found = 0
    sulfonamide_depth = float("inf")

    def dfs(node, depth=0):
        nonlocal sulfonamide_found, late_stage, max_depth_found, sulfonamide_depth, findings_json

        # Update max depth to determine late vs early stage
        max_depth_found = max(max_depth_found, depth)

        if node["type"] == "reaction":
            try:
                rsmi = node["metadata"]["mapped_reaction_smiles"]

                # Check for sulfonamide formation reactions
                for rxn in SULFONAMIDE_FORMATION_REACTIONS:
                    if checker.check_reaction(rxn, rsmi):
                        print(f"Found sulfonamide formation reaction at depth {depth}: {rsmi}")
                        sulfonamide_found = True
                        sulfonamide_depth = min(sulfonamide_depth, depth)
                        if rxn not in findings_json["atomic_checks"]["named_reactions"]:
                            findings_json["atomic_checks"]["named_reactions"].append(rxn)
            except Exception as e:
                print(f"Error in sulfonamide check: {e}")

        # Continue DFS traversal
        for child in node.get("children", []):
            # New logic for depth calculation:
            # Depth increases when going from chemical to reaction.
            # Depth stays the same when going from reaction to chemical.
            new_depth = depth
            if node["type"] != "reaction": # This means current node is 'chemical'
                new_depth = depth + 1
            
            dfs(child, new_depth)

    dfs(route)

    # Determine if sulfonamide formation is late-stage
    # Late stage is typically in the first third of the synthesis depth
    if max_depth_found > 0 and sulfonamide_found:
        current_late_stage = False
        if sulfonamide_depth <= 2:  # Absolute measure: depth <= 2 is late stage
            current_late_stage = True
        elif max_depth_found > 6:  # Relative measure for deeper trees
            current_late_stage = sulfonamide_depth <= max_depth_found // 3
        else:
            current_late_stage = sulfonamide_depth <= 2  # Default for medium-sized trees
        
        if current_late_stage:
            late_stage = True
            # Add structural constraint if late_stage condition is met
            findings_json["structural_constraints"].append({
                "type": "positional",
                "details": {
                    "targets": [
                        "Sulfonamide synthesis (Schotten-Baumann) primary amine",
                        "Sulfonamide synthesis (Schotten-Baumann) secondary amine"
                    ],
                    "position": "late_stage",
                    "condition": "The reaction must occur at a depth <= 2, or if the total route depth is greater than 6, at a depth less than or equal to one-third of the total depth."
                }
            })

    print(
        f"Late-stage sulfonamide formation detected: {late_stage} (sulfonamide found: {sulfonamide_found}, depth: {sulfonamide_depth}, max depth: {max_depth_found})"
    )
    return late_stage, findings_json
