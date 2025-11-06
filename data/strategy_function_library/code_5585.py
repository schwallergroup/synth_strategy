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
    This function detects late-stage alcohol protection via acetylation using acetic anhydride.
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

    found_acetylation = False
    acetylation_depth = None

    def is_alcohol_acetylation(rsmi):
        # The original function was severely flawed, containing logic for deprotection
        # and using incorrect checkers and string matching, leading to false positives.
        # This version is restricted to a single, high-confidence check to ensure accuracy.
        if checker.check_reaction("Acetic anhydride and alcohol to ester", rsmi):
            print(f"Found direct acetylation reaction")
            findings_json["atomic_checks"]["named_reactions"].append("Acetic anhydride and alcohol to ester")
            return True

        return False

    def dfs_traverse(node, depth=0):
        nonlocal found_acetylation, acetylation_depth, findings_json

        print(f"Traversing node at depth {depth}, type: {node['type']}")

        if node["type"] == "reaction":
            if "metadata" in node and "mapped_reaction_smiles" in node["metadata"]:
                rsmi = node["metadata"]["mapped_reaction_smiles"]

                if is_alcohol_acetylation(rsmi):
                    print(f"Found alcohol acetylation at depth {depth}")
                    found_acetylation = True
                    if acetylation_depth is None or depth < acetylation_depth:
                        acetylation_depth = depth

        # Traverse children
        for child in node.get("children", []):
            # New logic for depth calculation
            new_depth = depth
            if node['type'] != 'reaction': # If current node is chemical, depth increases
                new_depth = depth + 1
            # If current node is reaction, depth remains the same
            dfs_traverse(child, new_depth)

    # Start traversal
    dfs_traverse(route)

    print(
        f"Final result: found_acetylation={found_acetylation}, acetylation_depth={acetylation_depth}"
    )

    # Return True if acetylation was found at a late stage (depth <= 1)
    result = found_acetylation and acetylation_depth is not None and acetylation_depth <= 1

    if result:
        # Add the structural constraint if the overall condition is met
        findings_json["structural_constraints"].append({
            "type": "positional",
            "details": {
                "target": "Acetic anhydride and alcohol to ester",
                "constraint": {
                    "variable": "depth",
                    "operator": "<=",
                    "value": 1
                }
            }
        })

    return result, findings_json