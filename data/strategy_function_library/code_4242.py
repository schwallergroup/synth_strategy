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

def main(route) -> Tuple[bool, Dict]:
    """
    Detects the use of phthalimide as a protecting group for amines. This is identified by finding either a specific phthalimide deprotection step or a Mitsunobu reaction involving an imide, which is a common method for installing the phthalimide group (a variant of the Gabriel synthesis).
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

    found_phthalimide = False

    def dfs_traverse(node, depth=0):
        nonlocal found_phthalimide, findings_json

        if node["type"] == "reaction":
            if "mapped_reaction_smiles" in node.get("metadata", {}):
                try:
                    rsmi = node["metadata"]["mapped_reaction_smiles"]

                    print(f"Analyzing reaction: {rsmi}")

                    # Check if this is a phthalimide deprotection reaction
                    if checker.check_reaction("Phthalimide deprotection", rsmi):
                        print("Found phthalimide protection (via deprotection reaction)")
                        found_phthalimide = True
                        findings_json["atomic_checks"]["named_reactions"].append("Phthalimide deprotection")
                        # No return here, continue to check for Mitsunobu as well if needed, or if only one match is enough, then return

                    # Check for Mitsunobu reaction with phthalimide
                    if checker.check_reaction("Mitsunobu_imide", rsmi):
                        print("Found phthalimide protection via Mitsunobu reaction")
                        found_phthalimide = True
                        findings_json["atomic_checks"]["named_reactions"].append("Mitsunobu_imide")

                except Exception as e:
                    print(f"Error analyzing reaction: {e}")

        for child in node.get("children", []):
            # New logic for depth calculation
            new_depth = depth
            if node["type"] != "reaction": # If current node is chemical, depth increases
                new_depth = depth + 1
            # If current node is reaction, depth remains the same
            dfs_traverse(child, new_depth)

    dfs_traverse(route)

    # Add structural constraint if any of the target reactions were found
    if found_phthalimide:
        # This corresponds to the structural constraint in the metadata JSON
        # {"type": "count", "details": {"target": ["Phthalimide deprotection", "Mitsunobu_imide"], "operator": ">=", "value": 1}}
        findings_json["structural_constraints"].append({
            "type": "count",
            "details": {
                "target": [
                    "Phthalimide deprotection",
                    "Mitsunobu_imide"
                ],
                "operator": ">=",
                "value": 1
            }
        })

    print(f"Final result: {found_phthalimide}")
    return found_phthalimide, findings_json