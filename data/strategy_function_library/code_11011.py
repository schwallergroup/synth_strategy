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

# Combine for efficient checking
BOC_RELATED_REACTIONS = BOC_PROTECTION_REACTIONS + BOC_DEPROTECTION_REACTIONS

def main(route) -> Tuple[bool, Dict]:
    """
    Detects a Boc protecting group strategy (protection or deprotection) by matching the reaction to a predefined list of named reactions. The specific reactions checked are defined in the `BOC_RELATED_REACTIONS` list, which includes transformations like 'Boc amine protection' and 'Boc amine deprotection'.
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

    amine_protection_found = False
    boc_reaction_count = 0

    def dfs_traverse(node, depth=0):
        nonlocal amine_protection_found, boc_reaction_count, findings_json

        if node["type"] == "reaction" and "rsmi" in node.get("metadata", {}):
            rsmi = node["metadata"]["mapped_reaction_smiles"]
            print(f"Examining reaction: {rsmi}")

            # Check for any Boc-related protection or deprotection reactions
            for reaction_name in BOC_RELATED_REACTIONS:
                if checker.check_reaction(reaction_name, rsmi):
                    print(f"Detected {reaction_name} reaction")
                    amine_protection_found = True
                    boc_reaction_count += 1
                    if reaction_name not in findings_json["atomic_checks"]["named_reactions"]:
                        findings_json["atomic_checks"]["named_reactions"].append(reaction_name)
                    # No return here, continue to check for other reactions in case multiple apply

        for child in node.get("children", []):
            # New logic for depth calculation
            if node["type"] == "reaction":
                dfs_traverse(child, depth)
            else: # chemical node
                dfs_traverse(child, depth + 1)

    dfs_traverse(route)

    if boc_reaction_count >= 1:
        findings_json["structural_constraints"].append({
            "type": "count",
            "details": {
                "target": "Boc_protection_or_deprotection_reaction",
                "operator": ">=",
                "value": 1
            }
        })

    return amine_protection_found, findings_json
