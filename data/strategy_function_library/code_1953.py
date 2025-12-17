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
    Detects the use of silyl protection of hydroxyl groups
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

    silyl_protection_found = False

    def dfs_traverse(node, depth=0):
        nonlocal silyl_protection_found, findings_json

        if node["type"] == "reaction" and not silyl_protection_found:
            try:
                rsmi = node["metadata"]["mapped_reaction_smiles"]

                # Primary check: Use the reaction checker
                if checker.check_reaction("Alcohol protection with silyl ethers", rsmi):
                    silyl_protection_found = True
                    findings_json["atomic_checks"]["named_reactions"].append("Alcohol protection with silyl ethers")
                    print(f"Found silyl protection reaction: {rsmi}")
            except Exception as e:
                print(f"Error processing reaction node: {e}")

        # Traverse children
        for child in node.get("children", []):
            new_depth = depth
            if node["type"] != "reaction": # Chemical node
                new_depth = depth + 1
            # If node['type'] is 'reaction', new_depth remains 'depth'
            dfs_traverse(child, new_depth)

    # Start traversal
    dfs_traverse(route)

    if silyl_protection_found:
        findings_json["structural_constraints"].append({"type": "co-occurrence", "details": {"targets": ["Alcohol protection with silyl ethers"]}})

    return silyl_protection_found, findings_json
