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


# Refactored list of target reaction types
LATE_STAGE_ESTERIFICATION_TYPES = [
    "Esterification of Carboxylic Acids",
    "Schotten-Baumann to ester",
    "Transesterification",
    "O-alkylation of carboxylic acids with diazo compounds",
    "Oxidative esterification of primary alcohols",
]

def main(route) -> Tuple[bool, Dict]:
    """
    Detects if the synthesis route uses one of a specific list of named esterification reactions
    within the final two steps.
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

    late_stage_has_esterification = False

    def dfs_traverse(node, depth=0):
        nonlocal late_stage_has_esterification, findings_json

        # Consider final step and penultimate steps as "late-stage"
        if node["type"] == "reaction" and depth <= 2:
            if "rsmi" in node.get("metadata", {}):
                rsmi = node["metadata"]["mapped_reaction_smiles"]

                # Check for various esterification reactions
                for rxn_type in LATE_STAGE_ESTERIFICATION_TYPES:
                    if checker.check_reaction(rxn_type, rsmi):
                        late_stage_has_esterification = True
                        findings_json["atomic_checks"]["named_reactions"].append(rxn_type)
                        # If any esterification is found, we can add the structural constraint
                        # and then return, as the condition is 'any'.
                        findings_json["structural_constraints"].append({
                            "type": "positional",
                            "details": {
                                "targets": [
                                    "Esterification of Carboxylic Acids",
                                    "Schotten-Baumann to ester",
                                    "Transesterification",
                                    "O-alkylation of carboxylic acids with diazo compounds",
                                    "Oxidative esterification of primary alcohols"
                                ],
                                "position": "within_last_3_stages",
                                "condition": "any"
                            }
                        })
                        return

        for child in node.get("children", []):
            # New logic for depth calculation:
            # Depth increases only when traversing from a chemical node to a reaction node.
            # Depth remains the same when traversing from a reaction node to a chemical node.
            if node["type"] == "reaction":
                dfs_traverse(child, depth)
            else: # Assuming 'chemical' or other types that should increase depth
                dfs_traverse(child, depth + 1)

    dfs_traverse(route)
    return late_stage_has_esterification, findings_json
