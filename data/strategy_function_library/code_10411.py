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


ESTER_FORMATION_REACTIONS = [
    "Esterification of Carboxylic Acids",
    "O-alkylation of carboxylic acids with diazo compounds",
    "Oxidative esterification of primary alcohols",
    "Transesterification",
    "Schotten-Baumann to ester",
]

def main(route) -> Tuple[bool, Dict]:
    """
    Detects a synthetic strategy involving late-stage ester formation using one of the
    following reaction types: Esterification of Carboxylic Acids, O-alkylation of
    carboxylic acids with diazo compounds, Oxidative esterification of primary
    alcohols, Transesterification, Schotten-Baumann to ester.
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

    # Initialize tracking variable
    has_late_esterification = False

    def dfs_traverse(node, depth=0):
        nonlocal has_late_esterification, findings_json

        if node["type"] == "reaction":
            if "metadata" in node and "mapped_reaction_smiles" in node["metadata"]:
                rsmi = node["metadata"]["mapped_reaction_smiles"]

                try:
                    # Check for esterification reaction at late stage (depth <= 3)
                    if depth <= 3:
                        found_reaction_in_list = False
                        for r in ESTER_FORMATION_REACTIONS:
                            if checker.check_reaction(r, rsmi):
                                has_late_esterification = True
                                findings_json["atomic_checks"]["named_reactions"].append(r)
                                found_reaction_in_list = True
                                # No break here, as we want to record all found reactions from the list
                        
                        if found_reaction_in_list:
                            # Add the structural constraint if any of the target reactions are found at the specified depth
                            findings_json["structural_constraints"].append({
                                "type": "positional",
                                "details": {
                                    "target_group": [
                                        "Esterification of Carboxylic Acids",
                                        "O-alkylation of carboxylic acids with diazo compounds",
                                        "Oxidative esterification of primary alcohols",
                                        "Transesterification",
                                        "Schotten-Baumann to ester"
                                    ],
                                    "condition": "any",
                                    "position": {
                                        "variable": "depth",
                                        "operator": "<=",
                                        "value": 3
                                    }
                                }
                            })

                except Exception as e:
                    print(f"Error processing reaction: {e}")

        # Process children
        if not has_late_esterification:
            for child in node.get("children", []):
                # New logic for depth calculation
                if node["type"] == "reaction":
                    # Depth remains the same when traversing from reaction to chemical
                    dfs_traverse(child, depth)
                else:
                    # Depth increases when traversing from chemical to reaction
                    dfs_traverse(child, depth + 1)

    # Start traversal
    dfs_traverse(route)

    return has_late_esterification, findings_json