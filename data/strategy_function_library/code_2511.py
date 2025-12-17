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


SULFONAMIDE_REACTIONS_OF_INTEREST = [
    "Sulfonamide synthesis (Schotten-Baumann) primary amine",
    "Sulfonamide synthesis (Schotten-Baumann) secondary amine",
]

def main(route) -> Tuple[bool, Dict]:
    """
    Detects if the synthesis route uses late-stage sulfonamide formation via a predefined list of reaction types.
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

    sulfonamide_formation_depths = []
    result = False

    def dfs_traverse(node, current_depth=0):
        nonlocal result, findings_json
        if node["type"] == "reaction":
            if "mapped_reaction_smiles" in node.get("metadata", {}):
                rsmi = node["metadata"]["mapped_reaction_smiles"]

                # Check if this is a sulfonamide synthesis reaction using a predefined list of known reaction types
                is_sulfonamide_reaction = False
                for r_name in SULFONAMIDE_REACTIONS_OF_INTEREST:
                    if checker.check_reaction(r_name, rsmi):
                        is_sulfonamide_reaction = True
                        findings_json["atomic_checks"]["named_reactions"].append(r_name)
                        # No break here, as multiple reactions might match, and we want to record all of them

                if is_sulfonamide_reaction:
                    print(f"Found sulfonamide formation at depth: {current_depth}")
                    sulfonamide_formation_depths.append(current_depth)

        for child in node.get("children", []):
            # New logic for depth calculation
            if node["type"] == "reaction":
                # Depth remains the same when traversing from a reaction node to a chemical node
                dfs_traverse(child, current_depth)
            else:
                # Depth increases when traversing from a chemical node to a reaction node
                dfs_traverse(child, current_depth + 1)

    dfs_traverse(route)

    # Find the minimum depth where sulfonamide formation occurs
    min_depth = min(sulfonamide_formation_depths) if sulfonamide_formation_depths else None

    # Check if sulfonamide formation occurs at depth 0 or 1 (late stage)
    is_late_stage = min_depth is not None and min_depth <= 1
    print(f"Late-stage sulfonamide formation detected: {is_late_stage} (depth: {min_depth})")
    
    result = is_late_stage

    if result:
        # Add the structural constraint if late-stage sulfonamide formation is detected
        # This corresponds to the positional constraint in the strategy JSON
        findings_json["structural_constraints"].append({
            "type": "positional",
            "details": {
                "targets": [
                    "Sulfonamide synthesis (Schotten-Baumann) primary amine",
                    "Sulfonamide synthesis (Schotten-Baumann) secondary amine"
                ],
                "position": "late_stage (depth <= 1)"
            }
        })

    return result, findings_json
