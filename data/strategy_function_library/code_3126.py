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


# Refactoring for Enumeration: Isolate the lists of reaction types.
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
    Identifies a Boc protection strategy by checking for specific reaction types corresponding to the installation or removal of a Boc group. The function searches for reactions defined in the `BOC_PROTECTION_REACTIONS` and `BOC_DEPROTECTION_REACTIONS` lists.
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

    boc_protection_count = 0
    boc_deprotection_count = 0

    def dfs_traverse(node, depth=0):
        nonlocal boc_protection_count, boc_deprotection_count, findings_json

        if node["type"] == "reaction" and "metadata" in node and "mapped_reaction_smiles" in node["metadata"]:
            rsmi = node["metadata"]["mapped_reaction_smiles"]

            try:
                # Check for Boc protection reactions
                for reaction_type in BOC_PROTECTION_REACTIONS:
                    if checker.check_reaction(reaction_type, rsmi):
                        boc_protection_count += 1
                        findings_json["atomic_checks"]["named_reactions"].append(reaction_type)
                        break

                # Check for Boc deprotection reactions
                for reaction_type in BOC_DEPROTECTION_REACTIONS:
                    if checker.check_reaction(reaction_type, rsmi):
                        boc_deprotection_count += 1
                        findings_json["atomic_checks"]["named_reactions"].append(reaction_type)
                        break
            except Exception:
                # Silently ignore errors in reaction processing
                pass

        # Process children
        for child in node.get("children", []):
            # New logic for depth calculation
            new_depth = depth
            if node["type"] != "reaction": # If current node is chemical, depth increases
                new_depth = depth + 1
            # If current node is reaction, depth remains the same
            dfs_traverse(child, new_depth)

    # Start traversal
    dfs_traverse(route)

    # Determine the overall result
    result = boc_protection_count >= 1 or boc_deprotection_count >= 1

    # Add structural constraint if the overall condition is met
    if result:
        findings_json["structural_constraints"].append({
            "type": "count",
            "details": {
                "target": "Boc protection or deprotection",
                "operator": ">=",
                "value": 1
            }
        })

    # Return True if we have at least one Boc protection OR deprotection
    return result, findings_json
