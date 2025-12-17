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


LATE_STAGE_ETHERIFICATIONS = [
    "Williamson Ether Synthesis",
    "Mitsunobu aryl ether",
    "Chan-Lam etherification",
    "Williamson Ether Synthesis (intra to epoxy)",
    "Mitsunobu aryl ether (intramolecular)",
    "{Williamson ether}",
]

def main(route) -> Tuple[bool, Dict]:
    """
    Detects a convergent synthesis where fragments are joined via a late-stage etherification reaction from a predefined list of named reactions.
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

    # Track if we found the key features
    found_etherification = False
    fragment_branches = 0

    def dfs_traverse(node, depth=0):
        nonlocal found_etherification, fragment_branches, findings_json

        # Determine the depth for the next recursive call
        # Depth increases if current node is chemical (mol), stays same if reaction
        next_depth = depth + 1 if node["type"] == "mol" else depth

        # For reaction nodes, check for etherification at depth 0
        if node["type"] == "reaction":
            if depth == 0:
                rsmi = node.get("metadata", {}).get("rsmi", "")
                if rsmi:
                    # Check for known etherification reactions
                    for name in LATE_STAGE_ETHERIFICATIONS:
                        if checker.check_reaction(name, rsmi):
                            found_etherification = True
                            if name not in findings_json["atomic_checks"]["named_reactions"]:
                                findings_json["atomic_checks"]["named_reactions"].append(name)
                            print(f"Found late-stage etherification reaction at depth 0: {rsmi}")
                            break # Found one, no need to check others for this reaction node

            # Traverse children with appropriate depth
            for child in node.get("children", []):
                dfs_traverse(child, next_depth)

        # For molecule nodes
        else:  # node["type"] == "mol"
            # Count significant molecular fragments at depth 1
            # These are non-starting material molecules that represent key fragments
            if depth == 1 and not node.get("in_stock", False):
                fragment_branches += 1
                print(
                    f"Found fragment branch at depth 1: {node['smiles']}, total: {fragment_branches}"
                )

            # Traverse children with appropriate depth
            for child in node.get("children", []):
                dfs_traverse(child, next_depth)

    # Start traversal
    dfs_traverse(route)

    # Return True if we found both features of the strategy
    result = found_etherification and fragment_branches >= 2

    # Add structural constraints to findings_json if conditions are met
    if found_etherification:
        # This corresponds to the positional constraint for late-stage etherification
        findings_json["structural_constraints"].append({
            "type": "positional",
            "details": {
                "target": [
                    "Williamson Ether Synthesis",
                    "Mitsunobu aryl ether",
                    "Chan-Lam etherification",
                    "Williamson Ether Synthesis (intra to epoxy)",
                    "Mitsunobu aryl ether (intramolecular)",
                    "{Williamson ether}"
                ],
                "position": "last_stage"
            }
        })
    
    if fragment_branches >= 2:
        # This corresponds to the count constraint for non-stock precursor fragments
        findings_json["structural_constraints"].append({
            "type": "count",
            "details": {
                "target": "non-stock_precursor_fragments",
                "operator": ">=",
                "value": 2
            }
        })

    if found_etherification and fragment_branches >= 2:
        # This corresponds to the co-occurrence constraint
        findings_json["structural_constraints"].append({
            "type": "co-occurrence",
            "details": {
                "targets": [
                    "late-stage_etherification",
                    "convergent_fragment_count"
                ]
            }
        })

    print(f"Convergent synthesis with late-stage etherification: {result}")
    print(f"Found etherification: {found_etherification}, Fragment branches: {fragment_branches}")
    return result, findings_json
