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


SILYL_DEPROTECTION_REACTIONS = [
    "Alcohol deprotection from silyl ethers",
    "Alcohol deprotection from silyl ethers (double)",
    "Alcohol deprotection from silyl ethers (diol)",
]

def main(route) -> Tuple[bool, Dict]:
    """
    Detects a silyl ether protecting group strategy for alcohols. It identifies
    'Alcohol protection with silyl ethers' reactions and corresponding deprotection
    reactions, as defined in the SILYL_DEPROTECTION_REACTIONS list. A full
    strategy is flagged if a protection step is followed by a deprotection step
    later in the synthesis.
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

    result = False

    # Track protected molecules and their deprotection
    protected_molecules = {}  # product SMILES -> depth
    deprotected_molecules = {}  # product SMILES -> depth
    protection_reactions = []
    deprotection_reactions = []

    def dfs_traverse(node, depth=0):
        nonlocal result
        if node["type"] == "reaction" and "metadata" in node and "mapped_reaction_smiles" in node["metadata"]:
            rsmi = node["metadata"]["mapped_reaction_smiles"]
            product = rsmi.split(">")[-1]

            # Check for alcohol protection with silyl group
            if checker.check_reaction("Alcohol protection with silyl ethers", rsmi):
                protected_molecules[product] = depth
                protection_reactions.append(rsmi)
                if "Alcohol protection with silyl ethers" not in findings_json["atomic_checks"]["named_reactions"]:
                    findings_json["atomic_checks"]["named_reactions"].append("Alcohol protection with silyl ethers")

            # Check for alcohol deprotection from silyl ethers
            for name in SILYL_DEPROTECTION_REACTIONS:
                if checker.check_reaction(name, rsmi):
                    deprotected_molecules[product] = depth
                    deprotection_reactions.append(rsmi)
                    if name not in findings_json["atomic_checks"]["named_reactions"]:
                        findings_json["atomic_checks"]["named_reactions"].append(name)
                    break # Only need to find one deprotection reaction

        for child in node.get("children", []):
            # New logic for depth calculation
            new_depth = depth
            if node["type"] != "reaction": # If current node is chemical, depth increases
                new_depth = depth + 1
            # If current node is reaction, depth remains the same
            dfs_traverse(child, new_depth)

    dfs_traverse(route)

    # Check if we have both protection and deprotection steps
    if protection_reactions and deprotection_reactions:
        # In retrosynthesis, protection (earlier) should have higher depth than deprotection (later)
        for prot_depth in protected_molecules.values():
            for deprot_depth in deprotected_molecules.values():
                if (
                    prot_depth > deprot_depth
                ):  # Protection happened before deprotection in retrosynthesis
                    result = True
                    # Add structural constraint if both protection and deprotection are found in sequence
                    if {
                        "type": "sequence",
                        "details": {
                            "before": "Alcohol protection with silyl ethers",
                            "after": [
                                "Alcohol deprotection from silyl ethers",
                                "Alcohol deprotection from silyl ethers (double)",
                                "Alcohol deprotection from silyl ethers (diol)"
                            ]
                        }
                    } not in findings_json["structural_constraints"]:
                        findings_json["structural_constraints"].append({
                            "type": "sequence",
                            "details": {
                                "before": "Alcohol protection with silyl ethers",
                                "after": [
                                    "Alcohol deprotection from silyl ethers",
                                    "Alcohol deprotection from silyl ethers (double)",
                                    "Alcohol deprotection from silyl ethers (diol)"
                                ]
                            }
                        })
                    break # Found a valid sequence, no need to check further depths
            if result: # If result is True, break outer loop too
                break

    # If we have at least one protection step, consider it a partial strategy
    if protection_reactions and not result:
        result = True

    # If we have at least one deprotection step, it's still a silyl strategy
    if deprotection_reactions and not result:
        result = True

    return result, findings_json
