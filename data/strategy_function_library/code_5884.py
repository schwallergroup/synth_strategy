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


BOC_PROTECTION_REACTIONS = [
    "Boc amine protection",
    "Boc amine protection explicit",
    "Boc amine protection with Boc anhydride",
    "Boc amine protection (ethyl Boc)",
    "Boc amine protection of secondary amine",
    "Boc amine protection of primary amine",
]

BOC_DEPROTECTION_REACTIONS = [
    "Boc amine deprotection",
    "Boc amine deprotection of guanidine",
    "Boc amine deprotection to NH-NH2",
    "Tert-butyl deprotection of amine",
]

SULFONAMIDE_FORMATION_REACTIONS = [
    "Sulfonamide synthesis (Schotten-Baumann) primary amine",
    "Sulfonamide synthesis (Schotten-Baumann) secondary amine",
]

def main(route) -> Tuple[bool, Dict]:
    """
    Detects a multi-step synthetic strategy involving an optional Boc protection of an amine,
    followed by a Boc deprotection, and culminating in a late-stage (final or penultimate step)
    sulfonamide formation. The sequence is verified to occur in the correct order. The specific
    reaction types for each step are defined in the module-level constants:
    `BOC_PROTECTION_REACTIONS`, `BOC_DEPROTECTION_REACTIONS`, and `SULFONAMIDE_FORMATION_REACTIONS`.
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
    has_boc_protection = False
    has_boc_deprotection = False
    has_late_sulfonamide = False

    # Track depths for sequence verification
    boc_protection_depth = float("inf")
    boc_deprotection_depth = float("inf")
    sulfonamide_depth = float("inf")

    def dfs_traverse(node, depth=0):
        nonlocal has_boc_protection, has_boc_deprotection, has_late_sulfonamide
        nonlocal boc_protection_depth, boc_deprotection_depth, sulfonamide_depth
        nonlocal findings_json

        if node["type"] == "reaction":
            rsmi = node["metadata"]["mapped_reaction_smiles"]

            # Check for Boc protection
            for r in BOC_PROTECTION_REACTIONS:
                if checker.check_reaction(r, rsmi):
                    has_boc_protection = True
                    boc_protection_depth = min(boc_protection_depth, depth)
                    if r not in findings_json["atomic_checks"]["named_reactions"]:
                        findings_json["atomic_checks"]["named_reactions"].append(r)

            # Check for Boc deprotection
            for r in BOC_DEPROTECTION_REACTIONS:
                if checker.check_reaction(r, rsmi):
                    has_boc_deprotection = True
                    boc_deprotection_depth = min(boc_deprotection_depth, depth)
                    if r not in findings_json["atomic_checks"]["named_reactions"]:
                        findings_json["atomic_checks"]["named_reactions"].append(r)

            # Check for late-stage sulfonamide formation (depth 0 or 1)
            if depth <= 1:
                for r in SULFONAMIDE_FORMATION_REACTIONS:
                    if checker.check_reaction(r, rsmi):
                        has_late_sulfonamide = True
                        sulfonamide_depth = min(sulfonamide_depth, depth)
                        if r not in findings_json["atomic_checks"]["named_reactions"]:
                            findings_json["atomic_checks"]["named_reactions"].append(r)

        # Traverse children
        for child in node.get("children", []):
            # New logic for depth calculation
            new_depth = depth
            if node["type"] != "reaction": # If current node is chemical, depth increases
                new_depth = depth + 1
            # If current node is reaction, depth remains the same
            dfs_traverse(child, new_depth)

    # Start traversal from the root
    dfs_traverse(route)

    # Check if all conditions are met and in the correct sequence
    # We need at least deprotection and sulfonamide formation in the correct order
    # Protection is optional as it might be missing from the route
    correct_sequence = (
        has_boc_deprotection
        and has_late_sulfonamide
        and boc_deprotection_depth > sulfonamide_depth
    )

    # If we have protection, ensure it's in the correct sequence
    if has_boc_protection:
        correct_sequence = correct_sequence and (boc_deprotection_depth < boc_protection_depth)

    # Record structural constraints if conditions are met
    if has_late_sulfonamide and sulfonamide_depth <= 1:
        # Positional constraint for sulfonamide formation
        findings_json["structural_constraints"].append({
            "type": "positional",
            "details": {
                "target_reactions": [
                    "Sulfonamide synthesis (Schotten-Baumann) primary amine",
                    "Sulfonamide synthesis (Schotten-Baumann) secondary amine"
                ],
                "position": "late_stage",
                "max_depth": 1
            }
        })

    # Sequence constraint
    if has_boc_deprotection and has_late_sulfonamide and boc_deprotection_depth > sulfonamide_depth:
        sequence_constraint = {
            "type": "sequence",
            "details": {
                "event_definitions": {
                    "boc_protection": [
                        "Boc amine protection",
                        "Boc amine protection explicit",
                        "Boc amine protection with Boc anhydride",
                        "Boc amine protection (ethyl Boc)",
                        "Boc amine protection of secondary amine",
                        "Boc amine protection of primary amine"
                    ],
                    "boc_deprotection": [
                        "Boc amine deprotection",
                        "Boc amine deprotection of guanidine",
                        "Boc amine deprotection to NH-NH2",
                        "Tert-butyl deprotection of amine"
                    ],
                    "sulfonamide_formation": [
                        "Sulfonamide synthesis (Schotten-Baumann) primary amine",
                        "Sulfonamide synthesis (Schotten-Baumann) secondary amine"
                    ]
                },
                "ordered_events": [
                    {
                        "event": "boc_protection",
                        "optional": True
                    },
                    {
                        "event": "boc_deprotection",
                        "optional": False
                    },
                    {
                        "event": "sulfonamide_formation",
                        "optional": False
                    }
                ]
            }
        }
        # Only add if boc_protection is present or if it's optional and not required for the sequence
        # The logic for `correct_sequence` already handles the optionality of boc_protection
        if correct_sequence:
            findings_json["structural_constraints"].append(sequence_constraint)

    return correct_sequence, findings_json
