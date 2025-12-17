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


ESTERIFICATION_REACTIONS = [
    "Esterification of Carboxylic Acids",
    "O-alkylation of carboxylic acids with diazo compounds",
    "Transesterification",
    "Oxidative esterification of primary alcohols",
    "Schotten-Baumann to ester",
]

def main(route) -> Tuple[bool, Dict]:
    """
    Detects if the synthesis includes a late-stage esterification
    (conversion of carboxylic acid to ester in the final steps).
    """
    esterification_depth = None
    
    findings_template = {
        "atomic_checks": {
            "named_reactions": [],
            "ring_systems": [],
            "functional_groups": []
        },
        "structural_constraints": []
    }
    findings_json = copy.deepcopy(findings_template)

    def dfs_traverse(node, depth=0):
        nonlocal esterification_depth, findings_json

        if node["type"] == "reaction":
            if "metadata" in node and "mapped_reaction_smiles" in node["metadata"]:
                rsmi = node["metadata"]["mapped_reaction_smiles"]

                # Check for various esterification reactions
                for reaction_type in ESTERIFICATION_REACTIONS:
                    if checker.check_reaction(reaction_type, rsmi):
                        print(f"{reaction_type} detected at depth {depth}, rsmi: {rsmi}")
                        if esterification_depth is None or depth < esterification_depth:
                            esterification_depth = depth
                        findings_json["atomic_checks"]["named_reactions"].append(reaction_type)
                        break

        # Determine the new depth for the recursive call
        new_depth = depth
        if node["type"] != "reaction":  # If current node is 'chemical'
            new_depth = depth + 1

        # Continue traversing the tree
        for child in node.get("children", []):
            dfs_traverse(child, new_depth)

    # Start traversal from the root
    dfs_traverse(route)

    # Return True if esterification was found in the late stage (depth ≤ 3)
    result = esterification_depth is not None and esterification_depth <= 3
    print(f"Late stage esterification result: {result} (depth: {esterification_depth})")
    
    if result:
        # Add the structural constraint if the condition is met
        findings_json["structural_constraints"].append({
            "type": "positional",
            "details": {
                "targets": [
                    "Esterification of Carboxylic Acids",
                    "O-alkylation of carboxylic acids with diazo compounds",
                    "Transesterification",
                    "Oxidative esterification of primary alcohols",
                    "Schotten-Baumann to ester"
                ],
                "position": {
                    "operator": "<=",
                    "value": 3,
                    "unit": "depth_from_root"
                }
            }
        })

    return result, findings_json
