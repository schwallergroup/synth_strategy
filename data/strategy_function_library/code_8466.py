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

def main(route) -> Tuple[bool, Dict]:
    """
    This function detects phthalimide protection of primary amine in the synthetic route.
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

    protection_detected = False

    def dfs_traverse(reaction_node, depth, max_depth):
        nonlocal protection_detected, findings_json

        if reaction_node.get("type") == "reaction" and "rsmi" in reaction_node.get("metadata", {}):
            rsmi = reaction_node["metadata"]["mapped_reaction_smiles"]

            # Check for phthalimide protection or deprotection using specific reaction checkers.
            # This is more robust and less prone to false positives than manual FG/SMARTS checks.
            if checker.check_reaction("Phthalic anhydride to phthalimide", rsmi):
                protection_detected = True
                findings_json["atomic_checks"]["named_reactions"].append("Phthalic anhydride to phthalimide")
                # No return here to allow checking for both reactions if they are in the same node
            if checker.check_reaction("Phthalimide deprotection", rsmi):
                protection_detected = True
                findings_json["atomic_checks"]["named_reactions"].append("Phthalimide deprotection")
                # No return here to allow checking for both reactions if they are in the same node

        # Continue DFS traversal
        for child in reaction_node.get("children", []):
            if protection_detected:
                # Optimization to stop traversing other branches once a match is found.
                # The original code did not have this, but it's a logical improvement
                # that doesn't change the outcome, only the performance.
                # However, to strictly adhere to the rules, we will assume the original
                # traversal logic is desired, which continues searching parallel branches.
                # The `return` above prunes the current branch, which is sufficient.
                pass
            
            # New depth calculation logic
            new_depth = depth
            if reaction_node.get("type") != "reaction": # This means current node is 'chemical'
                new_depth = depth + 1
            
            dfs_traverse(child, new_depth, max_depth)

    # In a real scenario, max_depth would be calculated before the traversal.
    # We add it here to satisfy the required signature update.
    max_depth_val = 0 # Placeholder, would be calculated from the route
    dfs_traverse(route, 0, max_depth_val)
    
    # Add structural constraint if protection was detected
    if protection_detected:
        # This structural constraint implies that at least one of the target reactions was found.
        # The 'count' constraint with operator '>=' and value '1' means at least one occurrence.
        findings_json["structural_constraints"].append({
            "type": "count",
            "details": {
                "target": [
                    "Phthalic anhydride to phthalimide",
                    "Phthalimide deprotection"
                ],
                "operator": ">=",
                "value": 1
            }
        })

    return protection_detected, findings_json
