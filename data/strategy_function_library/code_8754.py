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


# Refactored lists of specific reaction types
PROTECTION_REACTIONS = [
    "Boc amine protection",
    "Boc amine protection explicit",
    "Boc amine protection with Boc anhydride",
    "Boc amine protection (ethyl Boc)",
    "Boc amine protection of secondary amine",
    "Boc amine protection of primary amine",
    "Alcohol protection with silyl ethers",
]

DEPROTECTION_REACTIONS = [
    "Phthalimide deprotection",
    "Boc amine deprotection",
    "Boc amine deprotection of guanidine",
    "Boc amine deprotection to NH-NH2",
    "Tert-butyl deprotection of amine",
    "Alcohol deprotection from silyl ethers",
    "Alcohol deprotection from silyl ethers (double)",
    "Alcohol deprotection from silyl ethers (diol)",
    "N-glutarimide deprotection",
    "Hydroxyl benzyl deprotection",
    "Carboxyl benzyl deprotection",
]

def main(route) -> Tuple[bool, Dict]:
    """
    Detects a potential protection-deprotection strategy by identifying specific named reactions from the PROTECTION_REACTIONS and DEPROTECTION_REACTIONS lists.
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
    protection_events = []
    deprotection_events = []

    def dfs_traverse(node, depth=0):
        nonlocal protection_events, deprotection_events, findings_json
        if node["type"] == "reaction" and "metadata" in node and "mapped_reaction_smiles" in node["metadata"]:
            rsmi = node["metadata"]["mapped_reaction_smiles"]

            # Check for protection reactions
            for reaction in PROTECTION_REACTIONS:
                if checker.check_reaction(reaction, rsmi):
                    protection_events.append((reaction, depth))
                    if reaction not in findings_json["atomic_checks"]["named_reactions"]:
                        findings_json["atomic_checks"]["named_reactions"].append(reaction)

            # Check for deprotection reactions
            for reaction in DEPROTECTION_REACTIONS:
                if checker.check_reaction(reaction, rsmi):
                    deprotection_events.append((reaction, depth))
                    if reaction not in findings_json["atomic_checks"]["named_reactions"]:
                        findings_json["atomic_checks"]["named_reactions"].append(reaction)

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

    result = False
    # Check if we have both protection and deprotection events
    if protection_events and deprotection_events:
        # Add co-occurrence constraint if met
        findings_json["structural_constraints"].append({
            "type": "co-occurrence",
            "details": {
                "targets": [
                    "any_protection_reaction",
                    "any_deprotection_reaction"
                ],
                "description": "The route must contain at least one reaction from the protection list and at least one reaction from the deprotection list."
            }
        })

        # Find the earliest protection and latest deprotection
        min_protection_depth = min(depth for _, depth in protection_events)
        max_deprotection_depth = max(depth for _, depth in deprotection_events)

        # In retrosynthetic traversal, protection should be at a lower depth than deprotection
        # (protection happens earlier in forward synthesis)
        valid_strategy = min_protection_depth < max_deprotection_depth
        if valid_strategy:
            findings_json["structural_constraints"].append({
                "type": "sequence",
                "details": {
                    "ordered_events": [
                        "any_protection_reaction",
                        "any_deprotection_reaction"
                    ],
                    "description": "The minimum retrosynthetic depth of a protection reaction must be less than the maximum retrosynthetic depth of a deprotection reaction, ensuring protection occurs earlier in the forward synthesis."
                }
            })
        result = valid_strategy

    return result, findings_json
