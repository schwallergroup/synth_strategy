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


ESTERIFICATION_REACTIONS = [
    "Esterification of Carboxylic Acids",
    "Transesterification",
    "Oxidative esterification of primary alcohols",
    "O-alkylation of carboxylic acids with diazo compounds",
]

def main(route) -> Tuple[bool, Dict]:
    """
    Checks for the repeated use of specific esterification reactions across a synthesis. This strategy is flagged if two or more reactions from the following list are found: Esterification of Carboxylic Acids, Transesterification, Oxidative esterification of primary alcohols, and O-alkylation of carboxylic acids with diazo compounds.
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

    # Count esterification reactions
    esterification_count = 0

    def dfs_traverse(node, depth=0):
        nonlocal esterification_count, findings_json

        if node["type"] == "reaction":
            if "rsmi" in node.get("metadata", {}):
                rsmi = node["metadata"]["mapped_reaction_smiles"]

                # Check if this is an esterification reaction using a predefined list of reaction types
                for reaction_name in ESTERIFICATION_REACTIONS:
                    if checker.check_reaction(reaction_name, rsmi):
                        esterification_count += 1
                        findings_json["atomic_checks"]["named_reactions"].append(reaction_name)
                        # Only count each reaction once per node, so break after first match
                        break 

        # Determine the next depth based on the current node's type
        next_depth = depth
        if node["type"] != "reaction": # If current node is chemical, depth increases
            next_depth = depth + 1

        # Process children
        for child in node.get("children", []):
            dfs_traverse(child, next_depth)

    # Start traversal
    dfs_traverse(route)

    # Strategy is present if multiple esterification reactions are found
    result = esterification_count >= 2

    if result:
        findings_json["structural_constraints"].append({
            "type": "count",
            "details": {
                "target": "esterification_reaction",
                "operator": ">=",
                "value": 2
            }
        })

    return result, findings_json
