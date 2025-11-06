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
    Checks if every intermediate molecule in the synthesis contains a carboxylic acid functional group. It does not check the final product or starting materials.
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

    carboxylic_acid_steps = []
    all_intermediates = []
    result = False

    def dfs_traverse(node, depth=0):
        nonlocal carboxylic_acid_steps, all_intermediates, findings_json

        if node["type"] == "mol":
            # Check if this is the final product (depth 0)
            if depth == 0:
                pass
            # Check if this is an intermediate molecule (not in stock)
            elif not node.get("in_stock", False):
                all_intermediates.append(depth)
                if checker.check_fg("Carboxylic acid", node["smiles"]):
                    carboxylic_acid_steps.append(depth)
                    findings_json["atomic_checks"]["functional_groups"].append("Carboxylic acid")

        # Determine the depth for the recursive call based on the current node's type
        new_depth = depth
        if node["type"] == "mol": # From chemical to reaction, depth increases
            new_depth = depth + 1
        # If node['type'] == 'reaction', depth remains the same (from reaction to chemical)

        # Traverse children
        for child in node.get("children", []):
            dfs_traverse(child, new_depth)

    # Start traversal from the root
    dfs_traverse(route)

    # Check if carboxylic acid is present in all intermediate steps
    if all_intermediates and carboxylic_acid_steps:
        # Check if all intermediate steps have carboxylic acid
        if set(all_intermediates) == set(carboxylic_acid_steps):
            result = True
            findings_json["structural_constraints"].append({
                "type": "positional",
                "details": {
                    "target": "Carboxylic acid",
                    "position": "all_intermediates"
                }
            })

    return result, findings_json
