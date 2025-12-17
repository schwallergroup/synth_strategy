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
    "Boc amine deprotection of guanidine",
    "Boc amine deprotection to NH-NH2",
    "Tert-butyl deprotection of amine",
]
SILYL_PROTECTION_REACTIONS = ["Alcohol protection with silyl ethers"]
SILYL_DEPROTECTION_REACTIONS = [
    "Alcohol deprotection from silyl ethers",
    "Alcohol deprotection from silyl ethers (double)",
    "Alcohol deprotection from silyl ethers (diol)",
]

ALL_PROTECTION_REACTIONS = BOC_PROTECTION_REACTIONS + SILYL_PROTECTION_REACTIONS
ALL_DEPROTECTION_REACTIONS = BOC_DEPROTECTION_REACTIONS + SILYL_DEPROTECTION_REACTIONS

def main(route) -> Tuple[bool, Dict]:
    """
    Detects if the synthetic route employs a protection-deprotection strategy for common protecting groups. It identifies protection and deprotection steps for Boc on amines and silyl ethers on alcohols, using a predefined list of specific reaction templates. A strategy is flagged if a deprotection step is found, or if a matched protection-deprotection pair is found in the correct synthetic order (protection before deprotection).
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

    protection_reactions = []
    deprotection_reactions = []
    result = False

    def dfs_traverse(node, depth=0):
        nonlocal result
        if node["type"] == "reaction":
            if "mapped_reaction_smiles" in node["metadata"]:
                rsmi = node["metadata"]["mapped_reaction_smiles"]

                # Check for protection reactions
                for r in ALL_PROTECTION_REACTIONS:
                    if checker.check_reaction(r, rsmi):
                        protection_reactions.append((rsmi, depth))
                        findings_json["atomic_checks"]["named_reactions"].append(r)
                        print(f"Protection reaction detected at depth {depth}: {rsmi}")

                # Check for deprotection reactions
                for r in ALL_DEPROTECTION_REACTIONS:
                    if checker.check_reaction(r, rsmi):
                        deprotection_reactions.append((rsmi, depth))
                        findings_json["atomic_checks"]["named_reactions"].append(r)
                        print(f"Deprotection reaction detected at depth {depth}: {rsmi}")

        # Continue traversing
        for child in node.get("children", []):
            # New logic for depth calculation
            new_depth = depth
            if node["type"] != "reaction": # If current node is chemical, depth increases
                new_depth = depth + 1
            # If current node is reaction, depth remains the same
            dfs_traverse(child, new_depth)

    # Start traversal
    dfs_traverse(route)

    print(f"Protection reactions found: {len(protection_reactions)}")
    print(f"Deprotection reactions found: {len(deprotection_reactions)}")

    # Check if we have both protection and deprotection
    if protection_reactions and deprotection_reactions:
        # In retrosynthetic traversal, protection happens at higher depth than deprotection
        for prot_rsmi, prot_depth in protection_reactions:
            for deprot_rsmi, deprot_depth in deprotection_reactions:
                # In retrosynthetic traversal, protection should be at higher depth
                if prot_depth > deprot_depth:
                    print(
                        f"Valid protection-deprotection sequence found: protection at depth {prot_depth}, deprotection at depth {deprot_depth}"
                    )
                    result = True
                    # Add structural constraint for sequence
                    findings_json["structural_constraints"].append({
                        "type": "sequence",
                        "details": {
                            "before": "any_protection_reaction",
                            "after": "any_deprotection_reaction",
                            "description": "A protection reaction from the list ALL_PROTECTION_REACTIONS must occur at a greater retrosynthetic depth (earlier in forward synthesis) than a deprotection reaction from ALL_DEPROTECTION_REACTIONS."
                        }
                    })
                    # If a valid sequence is found, we can break and return True
                    # However, the original logic continues to check the deprotection_reactions only case
                    # So, we set result and let the function proceed.

    # If we only have deprotection reactions, it's still a protection-deprotection strategy
    # since the protection might have happened outside the route
    if deprotection_reactions:
        print(
            "Only deprotection reactions found, but this still indicates a protection-deprotection strategy"
        )
        result = True
        # Add structural constraint for count
        findings_json["structural_constraints"].append({
            "type": "count",
            "details": {
                "target": "any_deprotection_reaction",
                "operator": ">=",
                "value": 1,
                "description": "The route is considered valid if at least one deprotection reaction from the list ALL_DEPROTECTION_REACTIONS is found, even without a corresponding protection step in the analyzed route."
            }
        })

    return result, findings_json
