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


AMIDE_FORMATION_REACTIONS = [
    "Acylation of Nitrogen Nucleophiles by Carboxylic Acids",
    "Acylation of Nitrogen Nucleophiles by Acyl/Thioacyl/Carbamoyl Halides and Analogs_N",
    "Acyl chloride with primary amine to amide (Schotten-Baumann)",
    "Acyl chloride with secondary amine to amide",
    "Carboxylic acid with primary amine to amide",
    "Ester with primary amine to amide",
    "Ester with secondary amine to amide",
    "Ester with ammonia to amide",
    "Acyl chloride with ammonia to amide",
    "Schotten-Baumann_amide",
]

def main(route) -> Tuple[bool, Dict]:
    """
    Detects a late-stage (final or penultimate step) amide formation strategy. The identification is based on matching the reaction to a predefined list of known amide formation reactions, such as acylation of amines with carboxylic acids, acyl halides, or esters.
    """
    late_amide_formation = False

    findings_template = {
        "atomic_checks": {
            "named_reactions": [],
            "ring_systems": [],
            "functional_groups": []
        },
        "structural_constraints": []
    }
    findings_json = copy.deepcopy(findings_template)

    # Define the full strategy JSON for reference to pick structural constraints
    full_strategy_json = {
      "function_id": "code_6617",
      "filepath": "../data/merged_good_perf/code_6617.py",
      "description": "Detects a late-stage (final or penultimate step) amide formation strategy. The identification is based on matching the reaction to a predefined list of known amide formation reactions, such as acylation of amines with carboxylic acids, acyl halides, or esters.",
      "atomic_checks": {
        "named_reactions": [
          "Acylation of Nitrogen Nucleophiles by Carboxylic Acids",
          "Acylation of Nitrogen Nucleophiles by Acyl/Thioacyl/Carbamoyl Halides and Analogs_N",
          "Acyl chloride with primary amine to amide (Schotten-Baumann)",
          "Acyl chloride with secondary amine to amide",
          "Carboxylic acid with primary amine to amide",
          "Ester with primary amine to amide",
          "Ester with secondary amine to amide",
          "Ester with ammonia to amide",
          "Acyl chloride with ammonia to amide",
          "Schotten-Baumann_amide"
        ],
        "ring_systems": [],
        "functional_groups": []
      },
      "structural_constraints": [
        {
          "type": "positional",
          "details": {
            "target": "amide_formation",
            "position": "last_or_penultimate_stage"
          }
        }
      ]
    }

    def dfs_traverse(node, depth=0):
        nonlocal late_amide_formation, findings_json

        if node["type"] == "reaction" and "metadata" in node and "mapped_reaction_smiles" in node["metadata"]:
            # Check for amide formation in late stage (depth <= 1)
            if depth <= 1:
                try:
                    rsmi = node["metadata"]["mapped_reaction_smiles"]
                    # Check if any known amide formation reaction occurred.
                    for rxn in AMIDE_FORMATION_REACTIONS:
                        if checker.check_reaction(rxn, rsmi):
                            late_amide_formation = True
                            findings_json["atomic_checks"]["named_reactions"].append(rxn)
                            # If a late-stage amide formation is found, add the structural constraint
                            for constraint in full_strategy_json["structural_constraints"]:
                                if constraint["type"] == "positional" and constraint["details"]["target"] == "amide_formation":
                                    findings_json["structural_constraints"].append(constraint)
                            return  # Strategy found, exit this traversal path
                except Exception:
                    # Silently ignore errors in reaction processing
                    pass

        # Recurse through children if the strategy has not been found yet.
        if not late_amide_formation:
            for child in node.get("children", []):
                # New logic for depth calculation
                if node["type"] == "reaction":
                    dfs_traverse(child, depth)
                else: # chemical node
                    dfs_traverse(child, depth + 1)

    dfs_traverse(route)
    return late_amide_formation, findings_json
