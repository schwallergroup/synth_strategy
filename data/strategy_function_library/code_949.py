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


PROTECTION_REACTIONS = {
    "silyl": ["Alcohol protection with silyl ethers", "Alcohol deprotection from silyl ethers"],
    "acetal/ketal": [
        "Aldehyde or ketone acetalization",
        "Acetal hydrolysis to aldehyde",
        "Ketal hydrolysis to ketone",
        "Acetal hydrolysis to diol",
    ],
    "boc": ["Boc amine protection", "Boc amine deprotection"],
    "benzyl": ["Hydroxyl benzyl deprotection", "Carboxyl benzyl deprotection"],
}

def main(route) -> Tuple[bool, Dict]:
    """
    Detects if a synthetic route employs two or more different protecting group strategies. A strategy is identified if any of its associated protection or deprotection reactions are found. The strategies checked are silyl, acetal/ketal, Boc, and benzyl, based on a predefined list of reaction names.
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

    # Dictionary to track protection strategies found
    protection_strategies = {
        "silyl": False,
        "acetal/ketal": False,  # For acetonide protection
        "boc": False,
        "benzyl": False,
    }

    def dfs_traverse(node, depth=0):
        nonlocal findings_json
        if node["type"] == "reaction" and "metadata" in node and "mapped_reaction_smiles" in node["metadata"]:
            rxn_smiles = node["metadata"]["mapped_reaction_smiles"]

            # Check for protection/deprotection reactions
            for strategy, reactions in PROTECTION_REACTIONS.items():
                for reaction in reactions:
                    if checker.check_reaction(reaction, rxn_smiles):
                        protection_strategies[strategy] = True
                        if reaction not in findings_json["atomic_checks"]["named_reactions"]:
                            findings_json["atomic_checks"]["named_reactions"].append(reaction)

        # Traverse children
        for child in node.get("children", []):
            # New logic for depth calculation
            if node["type"] == "chemical": # From chemical to reaction
                dfs_traverse(child, depth + 1)
            elif node["type"] == "reaction": # From reaction to chemical
                dfs_traverse(child, depth)
            else: # Default for other types, or if type is missing (shouldn't happen with current data)
                dfs_traverse(child, depth + 1)

    dfs_traverse(route)

    # Count how many different protection strategies were found
    strategies_found = sum(protection_strategies.values())

    # Determine the final result
    result = strategies_found >= 2

    # Add structural constraint if the condition is met
    if result:
        findings_json["structural_constraints"].append({
            "type": "count",
            "details": {
                "target": "protecting_group_strategies",
                "operator": ">=",
                "value": 2
            }
        })

    # Return True if at least 2 different protection strategies were found
    return result, findings_json
