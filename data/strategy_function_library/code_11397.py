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
    "Tert-butyl deprotection of amine",
]

def main(route) -> Tuple[bool, Dict]:
    """
    This function detects a synthetic strategy involving protection of an amine
    as a carbamate followed by transformation/deprotection.
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

    # Track protection and deprotection events
    protection_events = []
    deprotection_events = []

    result = False

    def dfs_traverse(node, depth=0):
        nonlocal result
        if node["type"] == "reaction":
            try:
                # Extract reactants and product from reaction SMILES
                rsmi = node["metadata"]["mapped_reaction_smiles"]

                # Check for Boc amine protection reactions
                is_boc_protection = False
                for r in BOC_PROTECTION_REACTIONS:
                    if checker.check_reaction(r, rsmi):
                        is_boc_protection = True
                        if r not in findings_json["atomic_checks"]["named_reactions"]:
                            findings_json["atomic_checks"]["named_reactions"].append(r)
                        break

                if is_boc_protection:
                    protection_events.append(depth)
                    print(f"Found amine protection at depth {depth}, reaction: {rsmi}")

                # Check for Boc deprotection
                is_boc_deprotection = False
                for r in BOC_DEPROTECTION_REACTIONS:
                    if checker.check_reaction(r, rsmi):
                        is_boc_deprotection = True
                        if r not in findings_json["atomic_checks"]["named_reactions"]:
                            findings_json["atomic_checks"]["named_reactions"].append(r)
                        break

                if is_boc_deprotection:
                    deprotection_events.append(depth)
                    print(f"Found carbamate deprotection at depth {depth}, reaction: {rsmi}")

            except Exception as e:
                print(f"Error processing reaction node: {e}")

        # Determine the next depth based on the current node's type
        next_depth = depth
        if node["type"] != "reaction":  # If current node is chemical, increase depth for children
            next_depth = depth + 1

        # Recursively process children
        for child in node.get("children", []):
            dfs_traverse(child, next_depth)

    # Start traversal from the root
    dfs_traverse(route)

    # Check if we have both protection and deprotection events in the correct order
    # (higher depth numbers come earlier in the synthesis)
    if protection_events and deprotection_events:
        print(f"Protection events at depths: {protection_events}")
        print(f"Deprotection events at depths: {deprotection_events}")
        for prot_depth in protection_events:
            for deprot_depth in deprotection_events:
                if (
                    prot_depth > deprot_depth
                ):  # Protection happens before deprotection in the synthesis
                    result = True
                    # Add the structural constraint to findings_json
                    structural_constraint_obj = {
                        "type": "sequence",
                        "details": {
                            "before": {
                                "event_type": "named_reaction",
                                "event_list": [
                                    "Boc amine protection",
                                    "Boc amine protection explicit",
                                    "Boc amine protection with Boc anhydride",
                                    "Boc amine protection (ethyl Boc)",
                                    "Boc amine protection of secondary amine",
                                    "Boc amine protection of primary amine"
                                ],
                                "description": "Protection of an amine with a Boc group."
                            },
                            "after": {
                                "event_type": "named_reaction",
                                "event_list": [
                                    "Boc amine deprotection",
                                    "Tert-butyl deprotection of amine"
                                ],
                                "description": "Deprotection of a Boc-protected amine."
                            }
                        }
                    }
                    if structural_constraint_obj not in findings_json["structural_constraints"]:
                        findings_json["structural_constraints"].append(structural_constraint_obj)
                    break # Found one valid sequence, no need to check other combinations
            if result: # If result is True, break outer loop as well
                break

    return result, findings_json
