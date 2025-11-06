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


LATE_STAGE_ESTERIFICATIONS = [
    "Esterification of Carboxylic Acids",
    "Transesterification",
    "Oxidative esterification of primary alcohols",
    "O-alkylation of carboxylic acids with diazo compounds",
    "Schotten-Baumann to ester",
]

def main(route) -> Tuple[bool, Dict]:
    """
    Detects if the synthesis route involves a late-stage esterification
    in the final step (depth=1), identified by checking if the reaction matches one of a specific list of named reactions.
    """
    result = False
    findings_template = {
        "atomic_checks": {
            "named_reactions": [],
            "ring_systems": [],
            "functional_groups": []
        },
        "structural_constraints": []
    }
    findings_json = copy.deepcopy(findings_template)

    # Define the full strategy JSON for structural constraints lookup
    strategy_json = {
      "function_id": "code_419",
      "filepath": "../data/merged_good_perf/code_419.py",
      "description": "Detects if the synthesis route involves a late-stage esterification in the final step (depth=1), identified by checking if the reaction matches one of a specific list of named reactions.",
      "atomic_checks": {
        "named_reactions": [
          "Esterification of Carboxylic Acids",
          "Transesterification",
          "Oxidative esterification of primary alcohols",
          "O-alkylation of carboxylic acids with diazo compounds",
          "Schotten-Baumann to ester"
        ],
        "ring_systems": [],
        "functional_groups": []
      },
      "structural_constraints": [
        {
          "type": "positional",
          "details": {
            "target": "Esterification of Carboxylic Acids",
            "position": "last_stage"
          }
        },
        {
          "type": "positional",
          "details": {
            "target": "Transesterification",
            "position": "last_stage"
          }
        },
        {
          "type": "positional",
          "details": {
            "target": "Oxidative esterification of primary alcohols",
            "position": "last_stage"
          }
        },
        {
          "type": "positional",
          "details": {
            "target": "O-alkylation of carboxylic acids with diazo compounds",
            "position": "last_stage"
          }
        },
        {
          "type": "positional",
          "details": {
            "target": "Schotten-Baumann to ester",
            "position": "last_stage"
          }
        }
      ]
    }

    def dfs_traverse(node, depth=0):
        nonlocal result, findings_json

        print(f"Examining node of type {node['type']} at depth {depth}")

        # Correctly check only the final reaction step (depth=1).
        if node["type"] == "reaction" and depth == 1:
            print(f"Examining potential late-stage reaction at depth {depth}")

            try:
                rsmi = node["metadata"]["mapped_reaction_smiles"]
                print(f"Reaction SMILES: {rsmi}")

                # Check for a list of specific, high-confidence esterification reactions.
                # This avoids false positives from simple FG checks.
                is_esterification = False
                for rxn in LATE_STAGE_ESTERIFICATIONS:
                    if checker.check_reaction(rxn, rsmi):
                        print(f"Found esterification reaction: {rsmi}")
                        is_esterification = True
                        result = True
                        findings_json["atomic_checks"]["named_reactions"].append(rxn)
                        # Add the corresponding structural constraint
                        for constraint in strategy_json["structural_constraints"]:
                            if constraint["details"]["target"] == rxn and constraint["details"]["position"] == "last_stage":
                                findings_json["structural_constraints"].append(constraint)
                        break # Found one, no need to check others

                if is_esterification:
                    return

            except Exception as e:
                print(f"Error analyzing reaction: {e}")

        # Continue traversing
        for child in node.get("children", []):
            # New logic for depth calculation:
            # Depth increases only when traversing from a chemical node to a reaction node.
            # Depth remains the same when traversing from a reaction node to a chemical node.
            new_depth = depth
            if node['type'] != 'reaction': # This means current node is 'chemical'
                new_depth = depth + 1
            dfs_traverse(child, new_depth)

    # Start traversal from the root
    dfs_traverse(route)
    print(f"Final result: {result}")
    return result, findings_json
