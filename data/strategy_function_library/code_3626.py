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


BOC_DEPROTECTION_REACTIONS = [
    "Boc amine deprotection",
    "Boc amine deprotection of guanidine",
    "Boc amine deprotection to NH-NH2",
    "Tert-butyl deprotection of amine",
]

def main(route) -> Tuple[bool, Dict]:
    """
    This function detects if the synthetic route involves a specific Boc deprotection reaction in the final step (depth=1).
    It checks for the following reaction types: 'Boc amine deprotection', 'Boc amine deprotection of guanidine', 'Boc amine deprotection to NH-NH2', and 'Tert-butyl deprotection of amine'.
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

    min_depth_deprotection = float("inf")
    has_deprotection = False

    def dfs_traverse(node, depth=0):
        nonlocal min_depth_deprotection, has_deprotection, findings_json

        # Process current node
        if node["type"] == "reaction":
            rsmi = node["metadata"].get("rsmi", "")

            # Check if this is a Boc deprotection reaction
            is_boc_deprotection = False

            for rxn_type in BOC_DEPROTECTION_REACTIONS:
                if checker.check_reaction(rxn_type, rsmi):
                    is_boc_deprotection = True
                    has_deprotection = True
                    min_depth_deprotection = min(min_depth_deprotection, depth)
                    findings_json["atomic_checks"]["named_reactions"].append(rxn_type)
                    print(f"Found {rxn_type} at depth {depth}")
                    break

        # Traverse children
        for child in node.get("children", []):
            # New logic for depth calculation
            new_depth = depth
            if node["type"] != "reaction": # If current node is chemical, depth increases
                new_depth = depth + 1
            # If current node is reaction, depth remains the same
            dfs_traverse(child, new_depth)

    # Start traversal from the root
    dfs_traverse(route)

    # Check if deprotection occurs at the latest stage (depth 1)
    is_late_stage = has_deprotection and min_depth_deprotection == 1

    if is_late_stage:
        findings_json["structural_constraints"].append({
            "type": "positional",
            "details": {
                "targets": [
                    "Boc amine deprotection",
                    "Boc amine deprotection of guanidine",
                    "Boc amine deprotection to NH-NH2",
                    "Tert-butyl deprotection of amine"
                ],
                "position": "last_stage"
            }
        })

    print(f"Has Boc deprotection: {has_deprotection}")
    print(f"Minimum depth of deprotection: {min_depth_deprotection}")
    print(f"Is late stage deprotection: {is_late_stage}")

    return is_late_stage, findings_json
