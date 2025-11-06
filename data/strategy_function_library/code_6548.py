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


CS_FORMING_REACTIONS = [
    "S-alkylation of thiols",
    "S-alkylation of thiols (ethyl)",
    "S-alkylation of thiols with alcohols",
    "S-alkylation of thiols with alcohols (ethyl)",
    "thia-Michael addition",
    "Newman-Kwart rearrangement",
]

def main(route) -> Tuple[bool, Dict]:
    """
    Checks if the synthesis route contains at least one C-S bond formation reaction from a predefined list: S-alkylation of thiols, S-alkylation of thiols (ethyl), S-alkylation of thiols with alcohols, S-alkylation of thiols with alcohols (ethyl), thia-Michael addition, or Newman-Kwart rearrangement.
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

    cs_bond_formation_count = 0

    def dfs_traverse(node, depth=0):
        nonlocal cs_bond_formation_count, findings_json

        if node["type"] == "reaction":
            if "rsmi" in node.get("metadata", {}):
                rsmi = node["metadata"]["mapped_reaction_smiles"]

                # Check for specific C-S bond forming reactions
                for reaction_type in CS_FORMING_REACTIONS:
                    if checker.check_reaction(reaction_type, rsmi):
                        print(
                            f"Detected C-S bond formation reaction '{reaction_type}' at depth {depth}"
                        )
                        cs_bond_formation_count += 1
                        findings_json["atomic_checks"]["named_reactions"].append(reaction_type)
                        break

        # Process children
        for child in node.get("children", []):
            # New logic for depth calculation
            if node["type"] == "reaction":
                dfs_traverse(child, depth)
            else:  # node['type'] == 'chemical'
                dfs_traverse(child, depth + 1)

    # Start traversal
    dfs_traverse(route)
    print(f"Total C-S bond formations detected: {cs_bond_formation_count}")

    result = cs_bond_formation_count >= 1

    if result:
        findings_json["structural_constraints"].append({
            "type": "count",
            "details": {
                "target": "C-S bond formation",
                "operator": ">=",
                "value": 1
            }
        })

    return result, findings_json
