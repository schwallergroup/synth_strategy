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

def main(route) -> Tuple[bool, Dict]:
    """
    This function detects if the synthesis route involves a Boc protection/deprotection sequence.
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

    # Track protection and deprotection events with their depths
    protection_events = []
    deprotection_events = []

    def dfs_traverse(node, depth=0):
        nonlocal protection_events, deprotection_events, findings_json
        if node["type"] == "reaction" and "metadata" in node and "mapped_reaction_smiles" in node["metadata"]:
            rsmi = node["metadata"]["mapped_reaction_smiles"]
            print(f"Depth {depth}, Examining reaction: {rsmi}")

            # Check for Boc protection reactions
            for r in BOC_PROTECTION_REACTIONS:
                if checker.check_reaction(r, rsmi):
                    print(f"Depth {depth}: Detected Boc protection reaction")
                    protection_events.append((depth, rsmi))
                    if r not in findings_json["atomic_checks"]["named_reactions"]:
                        findings_json["atomic_checks"]["named_reactions"].append(r)
                    break # Only need to find one match

            # Check for Boc deprotection reactions
            for r in BOC_DEPROTECTION_REACTIONS:
                if checker.check_reaction(r, rsmi):
                    print(f"Depth {depth}: Detected Boc deprotection reaction")
                    deprotection_events.append((depth, rsmi))
                    if r not in findings_json["atomic_checks"]["named_reactions"]:
                        findings_json["atomic_checks"]["named_reactions"].append(r)
                    break # Only need to find one match

        # Determine the next depth based on the current node's type
        next_depth = depth
        if node["type"] != "reaction": # This means it's a chemical node
            next_depth = depth + 1

        # Traverse children
        for child in node.get("children", []):
            dfs_traverse(child, next_depth)

    # Start traversal
    dfs_traverse(route)

    # Check if both protection and deprotection occurred
    has_protection = len(protection_events) > 0
    has_deprotection = len(deprotection_events) > 0

    print(f"Boc protection events: {len(protection_events)}")
    for depth, rsmi in protection_events:
        print(f"  Depth {depth}: {rsmi}")

    print(f"Boc deprotection events: {len(deprotection_events)}")
    for depth, rsmi in deprotection_events:
        print(f"  Depth {depth}: {rsmi}")

    # In retrosynthetic analysis, deprotection should be encountered at a lower depth than protection
    # (since we're traversing backwards from the target)
    has_correct_sequence = False
    if has_protection and has_deprotection:
        # Add the co-occurrence constraint if both are found
        findings_json["structural_constraints"].append({
            "type": "co-occurrence",
            "details": {
                "targets": [
                    "any_Boc_protection_reaction",
                    "any_Boc_deprotection_reaction"
                ]
            }
        })

        # Get the minimum depth for each event type
        min_protection_depth = min([depth for depth, _ in protection_events])
        min_deprotection_depth = min([depth for depth, _ in deprotection_events])

        # In retrosynthesis, deprotection should be encountered before protection
        # (lower depth = later stage in forward synthesis)
        if min_deprotection_depth < min_protection_depth:
            has_correct_sequence = True
            print(
                f"Correct sequence: Deprotection (depth {min_deprotection_depth}) occurs before Protection (depth {min_protection_depth}) in retrosynthetic analysis"
            )
            # Add the sequence constraint if the order is correct
            findings_json["structural_constraints"].append({
                "type": "sequence",
                "details": {
                    "before": "any_Boc_deprotection_reaction",
                    "after": "any_Boc_protection_reaction"
                }
            })
        else:
            print(
                f"Incorrect sequence: Protection (depth {min_protection_depth}) occurs before Deprotection (depth {min_deprotection_depth}) in retrosynthetic analysis"
            )

    result = has_protection and has_deprotection and has_correct_sequence
    print(f"Boc protection/deprotection sequence detected: {result}")
    return result, findings_json
