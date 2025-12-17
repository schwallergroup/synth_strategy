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
    Detects if a silyl protecting group strategy for alcohols is used in the synthesis.
    It identifies a relevant step if a reaction is classified as 'Alcohol protection with silyl ethers'
    or any of the specified deprotection reactions in SILYL_DEPROTECTION_REACTIONS.
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

    found_silyl_protection = False

    def dfs_traverse(node, depth=0):
        nonlocal found_silyl_protection, findings_json

        if node["type"] == "reaction" and not found_silyl_protection:
            rsmi = node["metadata"].get("rsmi")
            if not rsmi:
                return

            # Check for explicit silyl protection reaction
            if checker.check_reaction("Alcohol protection with silyl ethers", rsmi):
                found_silyl_protection = True
                findings_json["atomic_checks"]["named_reactions"].append("Alcohol protection with silyl ethers")
                return

            # Check for deprotection reactions
            for rxn in SILYL_DEPROTECTION_REACTIONS:
                if checker.check_reaction(rxn, rsmi):
                    found_silyl_protection = True
                    findings_json["atomic_checks"]["named_reactions"].append(rxn)
                    return

        # Recursively process children
        for child in node.get("children", []):
            if not found_silyl_protection:
                # New logic for depth calculation
                if node["type"] == "reaction":
                    dfs_traverse(child, depth)
                else: # node['type'] == 'chemical'
                    dfs_traverse(child, depth + 1)

    dfs_traverse(route)

    if found_silyl_protection:
        # This structural constraint is met if any relevant reaction is found
        findings_json["structural_constraints"].append({
            "type": "co-occurrence",
            "details": {
                "targets": [
                    "Alcohol protection with silyl ethers",
                    "Alcohol deprotection from silyl ethers",
                    "Alcohol deprotection from silyl ethers (double)",
                    "Alcohol deprotection from silyl ethers (diol)"
                ],
                "min_occurrences": 1,
                "description": "The synthesis route must contain at least one of the specified silyl protection or deprotection reactions."
            }
        })

    return found_silyl_protection, findings_json
