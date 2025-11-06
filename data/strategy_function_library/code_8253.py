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
    This function detects if the synthesis uses a ketal protection/deprotection strategy
    for a carbonyl group.
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

    found_ketal_formation = False
    found_ketal_deprotection = False

    def dfs_traverse(node, depth=0):
        nonlocal found_ketal_formation, found_ketal_deprotection, findings_json

        if node["type"] == "reaction" and "metadata" in node and "mapped_reaction_smiles" in node["metadata"]:
            rsmi = node["metadata"]["mapped_reaction_smiles"]

            # Check for acetal/ketal formation (protection)
            if checker.check_reaction("Aldehyde or ketone acetalization", rsmi):
                print(f"Found ketal formation at depth {depth}, rsmi: {rsmi}")
                found_ketal_formation = True
                findings_json["atomic_checks"]["named_reactions"].append("Aldehyde or ketone acetalization")

            # Check for acetal/ketal deprotection
            if checker.check_reaction("Ketal hydrolysis to ketone", rsmi):
                print(f"Found ketal deprotection at depth {depth}, rsmi: {rsmi}")
                found_ketal_deprotection = True
                findings_json["atomic_checks"]["named_reactions"].append("Ketal hydrolysis to ketone")

        # Traverse children
        for child in node.get("children", []):
            # New logic: depth increases only when traversing from chemical to reaction
            new_depth = depth
            if node["type"] != "reaction": # This means current node is 'chemical'
                new_depth = depth + 1
            dfs_traverse(child, new_depth)

    # Start traversal from root
    dfs_traverse(route)

    # Return True only if both protection and deprotection are found
    result = found_ketal_formation and found_ketal_deprotection
    print(f"Ketal protection strategy detected: {result}")
    print(f"Found ketal formation: {found_ketal_formation}")
    print(f"Found ketal deprotection: {found_ketal_deprotection}")

    # Add structural constraint if both conditions are met
    if result:
        findings_json["structural_constraints"].append({
            "type": "co-occurrence",
            "details": {
                "targets": [
                    "Aldehyde or ketone acetalization",
                    "Ketal hydrolysis to ketone"
                ]
            }
        })

    return result, findings_json
