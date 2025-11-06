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


BORYLATION_REACTION_TYPES = [
    "Preparation of boronic acids",
    "Preparation of boronic esters"
]

def main(route) -> Tuple[bool, Dict]:
    """
    This function detects if the synthesis route uses borylation of aryl halides
    to prepare coupling partners for cross-coupling reactions.
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

    borylation_found = False

    def dfs_traverse(node, depth=0):
        nonlocal borylation_found, findings_json

        if node["type"] == "reaction" and "metadata" in node and "mapped_reaction_smiles" in node["metadata"]:
            rsmi = node["metadata"]["mapped_reaction_smiles"]

            # Check if the reaction is one of the specified borylation types
            for r_type in BORYLATION_REACTION_TYPES:
                if checker.check_reaction(r_type, rsmi):
                    borylation_found = True
                    if r_type not in findings_json["atomic_checks"]["named_reactions"]:
                        findings_json["atomic_checks"]["named_reactions"].append(r_type)
                    break

        for child in node.get("children", []):
            if borylation_found:
                break
            
            # New depth calculation logic
            new_depth = depth
            if node["type"] != "reaction": # chemical node
                new_depth = depth + 1
            
            dfs_traverse(child, new_depth)

    dfs_traverse(route)
    
    if borylation_found:
        # Add the structural constraint if borylation was found
        findings_json["structural_constraints"].append({
            "type": "count",
            "details": {
                "target": [
                    "Preparation of boronic acids",
                    "Preparation of boronic esters"
                ],
                "operator": ">=",
                "value": 1
            }
        })

    return borylation_found, findings_json
