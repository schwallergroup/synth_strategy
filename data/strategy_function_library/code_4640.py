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
    Detects the use of Boc protection or deprotection of amines by checking for a specific list of named reactions. The function iterates through all reaction steps and returns True if any reaction matches a name in the predefined lists of protection or deprotection reactions.
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

    has_boc_protection = False

    def dfs_traverse(node, depth=0):
        nonlocal has_boc_protection, findings_json

        if has_boc_protection:
            return  # Early exit if we already found Boc protection

        if node["type"] == "reaction":
            # Extract reactants and product
            rsmi = node["metadata"]["mapped_reaction_smiles"]
            print(f"Depth {depth}, Examining reaction: {rsmi}")

            # Check for protection reactions
            for reaction_type in BOC_PROTECTION_REACTIONS:
                if checker.check_reaction(reaction_type, rsmi):
                    print(f"Found {reaction_type} reaction")
                    has_boc_protection = True
                    findings_json["atomic_checks"]["named_reactions"].append(reaction_type)
                    # Add structural constraint if found
                    findings_json["structural_constraints"].append({
                        "type": "count",
                        "details": {
                            "target": "Boc protection or deprotection reaction",
                            "operator": ">=",
                            "value": 1
                        }
                    })
                    return

            # Check for deprotection reactions
            for reaction_type in BOC_DEPROTECTION_REACTIONS:
                if checker.check_reaction(reaction_type, rsmi):
                    print(f"Found {reaction_type} reaction")
                    has_boc_protection = True
                    findings_json["atomic_checks"]["named_reactions"].append(reaction_type)
                    # Add structural constraint if found
                    findings_json["structural_constraints"].append({
                        "type": "count",
                        "details": {
                            "target": "Boc protection or deprotection reaction",
                            "operator": ">=",
                            "value": 1
                        }
                    })
                    return

        # Determine the next depth based on the current node's type
        next_depth = depth
        if node["type"] != "reaction": # If current node is chemical, depth increases
            next_depth = depth + 1

        # Traverse children
        for child in node.get("children", []):
            dfs_traverse(child, next_depth)

    # Start traversal from root
    dfs_traverse(route)

    print(f"Final result: has_boc_protection = {has_boc_protection}")
    return has_boc_protection, findings_json
