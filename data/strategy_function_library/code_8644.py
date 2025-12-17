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


# Refactoring for Enumeration: Isolate the lists of reaction names
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
    Detects the use of a Boc protection strategy by identifying specific Boc protection or deprotection reactions. The reactions checked are defined in the `BOC_PROTECTION_REACTIONS` and `BOC_DEPROTECTION_REACTIONS` lists.
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

        # Check for Boc protection/deprotection reactions
        if node["type"] == "reaction":
            try:
                rsmi = node["metadata"]["mapped_reaction_smiles"]

                # Check for Boc protection reactions
                for name in BOC_PROTECTION_REACTIONS:
                    if checker.check_reaction(name, rsmi):
                        has_boc_protection = True
                        if name not in findings_json["atomic_checks"]["named_reactions"]:
                            findings_json["atomic_checks"]["named_reactions"].append(name)
                        print(f"Detected Boc protection reaction at depth {depth}")
                        # No break here, as multiple reactions might match, and we want to record all

                # Check for Boc deprotection reactions
                for name in BOC_DEPROTECTION_REACTIONS:
                    if checker.check_reaction(name, rsmi):
                        has_boc_protection = True
                        if name not in findings_json["atomic_checks"]["named_reactions"]:
                            findings_json["atomic_checks"]["named_reactions"].append(name)
                        print(f"Detected Boc deprotection reaction at depth {depth}")
                        # No break here, as multiple reactions might match, and we want to record all

            except Exception as e:
                print(f"Error processing reaction node: {e}")

        # Determine the next depth based on the current node's type
        next_depth = depth
        if node["type"] != "reaction": # Depth increases only when traversing from a chemical node to a reaction node
            next_depth = depth + 1

        # Traverse children
        for child in node.get("children", []):
            dfs_traverse(child, next_depth)

    # Start traversal
    dfs_traverse(route)

    if has_boc_protection:
        # Add the structural constraint if any Boc protection/deprotection reaction was found
        # This corresponds to the strategy's structural_constraints entry:
        # {"type": "count", "details": {"target": "boc_protection_or_deprotection_reaction", "operator": ">=", "value": 1}}
        findings_json["structural_constraints"].append({
            "type": "count",
            "details": {
                "target": "boc_protection_or_deprotection_reaction",
                "operator": ">=",
                "value": 1
            }
        })

    print(f"Boc protection strategy: {has_boc_protection}")
    return has_boc_protection, findings_json
