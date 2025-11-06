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
    This function detects if a Boc protecting group is present in the final product and in the products of at least two other reaction steps in the synthesis.
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

    boc_present = []

    def dfs_traverse(node, depth=0):
        nonlocal boc_present, findings_json

        if node["type"] == "reaction":
            # Extract product
            rsmi = node["metadata"]["mapped_reaction_smiles"]
            product = rsmi.split(">")[-1]

            # Check for Boc group in product
            if checker.check_fg("Boc", product):
                boc_present.append(depth)
                if "Boc" not in findings_json["atomic_checks"]["functional_groups"]:
                    findings_json["atomic_checks"]["functional_groups"].append("Boc")
                print(f"Found Boc group in reaction product at depth {depth}")

        # Traverse children
        for child in node.get("children", []):
            new_depth = depth
            if node["type"] != "reaction": # If current node is chemical, depth increases
                new_depth = depth + 1
            # If current node is reaction, depth remains the same
            dfs_traverse(child, new_depth)

    # Call dfs_traverse on the root node
    dfs_traverse(route)

    # Check if Boc is present throughout synthesis (at the final product and at least at 3 different depths)
    strategy_detected = 0 in boc_present and len(boc_present) >= 3
    print(f"Boc-protected amine strategy detected: {strategy_detected}")
    print(f"Boc found at depths: {sorted(boc_present)}")

    if strategy_detected:
        if 0 in boc_present:
            findings_json["structural_constraints"].append({
                "type": "positional",
                "details": {
                    "target": "Boc",
                    "position": "last_stage"
                }
            })
        if len(boc_present) >= 3:
            findings_json["structural_constraints"].append({
                "type": "count",
                "details": {
                    "target": "Boc",
                    "operator": ">=",
                    "value": 3
                }
            })

    return strategy_detected, findings_json
