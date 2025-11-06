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


HYDROLYSIS_REACTIONS_OF_INTEREST = [
    "Hydrolysis or Hydrogenolysis of Carboxylic Esters or Thioesters",
    "Ester saponification (methyl deprotection)",
    "Ester saponification (alkyl deprotection)",
    "COOH ethyl deprotection",
]
ESTERIFICATION_REACTIONS_OF_INTEREST = [
    "Esterification of Carboxylic Acids",
    "O-alkylation of carboxylic acids with diazo compounds",
    "Transesterification",
]

def main(route) -> Tuple[bool, Dict]:
    """
    Checks for late-stage esterification or its retrosynthetic equivalent (hydrolysis)
    in the final or penultimate step. This is determined by checking for specific
    reaction types, including those in HYDROLYSIS_REACTIONS_OF_INTEREST and
    ESTERIFICATION_REACTIONS_OF_INTEREST.
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

    esterification_detected = False

    # Define the full strategy JSON for structural constraint lookup
    full_strategy_json = {
      "function_id": "code_10296",
      "filepath": "../data/merged_good_perf/code_10296.py",
      "description": "Checks for late-stage esterification or its retrosynthetic equivalent (hydrolysis) within the first three retrosynthetic steps. This is determined by checking for specific named hydrolysis or esterification reactions that are accompanied by the expected functional group transformation (ester to acid, or acid to ester).",
      "atomic_checks": {
        "named_reactions": [
          "Hydrolysis or Hydrogenolysis of Carboxylic Esters or Thioesters",
          "Ester saponification (methyl deprotection)",
          "Ester saponification (alkyl deprotection)",
          "COOH ethyl deprotection",
          "Esterification of Carboxylic Acids",
          "O-alkylation of carboxylic acids with diazo compounds",
          "Transesterification"
        ],
        "ring_systems": [],
        "functional_groups": [
          "Ester",
          "Carboxylic acid"
        ]
      },
      "structural_constraints": [
        {
          "type": "positional",
          "details": {
            "target": "reaction_evaluation",
            "position": "within_first_3_retrosynthetic_steps"
          }
        },
        {
          "type": "co-occurrence",
          "details": {
            "targets": [
              "Hydrolysis or Hydrogenolysis of Carboxylic Esters or Thioesters",
              "Ester",
              "Carboxylic acid"
            ]
          }
        },
        {
          "type": "co-occurrence",
          "details": {
            "targets": [
              "Ester saponification (methyl deprotection)",
              "Ester",
              "Carboxylic acid"
            ]
          }
        },
        {
          "type": "co-occurrence",
          "details": {
            "targets": [
              "Ester saponification (alkyl deprotection)",
              "Ester",
              "Carboxylic acid"
            ]
          }
        },
        {
          "type": "co-occurrence",
          "details": {
            "targets": [
              "COOH ethyl deprotection",
              "Ester",
              "Carboxylic acid"
            ]
          }
        },
        {
          "type": "co-occurrence",
          "details": {
            "targets": [
              "Esterification of Carboxylic Acids",
              "Carboxylic acid",
              "Ester"
            ]
          }
        },
        {
          "type": "co-occurrence",
          "details": {
            "targets": [
              "O-alkylation of carboxylic acids with diazo compounds",
              "Carboxylic acid",
              "Ester"
            ]
          }
        },
        {
          "type": "co-occurrence",
          "details": {
            "targets": [
              "Transesterification",
              "Carboxylic acid",
              "Ester"
            ]
          }
        }
      ]
    }

    def dfs_traverse(node, depth=0):
        nonlocal esterification_detected, findings_json

        if node["type"] == "reaction" and depth <= 2:  # Final or penultimate step
            # Record positional constraint if met
            if depth <= 2:
                positional_constraint = next((c for c in full_strategy_json["structural_constraints"] if c.get("type") == "positional" and c["details"].get("position") == "within_first_3_retrosynthetic_steps"), None)
                if positional_constraint and positional_constraint not in findings_json["structural_constraints"]:
                    findings_json["structural_constraints"].append(copy.deepcopy(positional_constraint))

            try:
                rsmi = node["metadata"]["mapped_reaction_smiles"]
                reactants = rsmi.split(">")[0].split(".")
                product = rsmi.split(">")[-1]

                is_hydrolysis = False
                for reaction_name in HYDROLYSIS_REACTIONS_OF_INTEREST:
                    if checker.check_reaction(reaction_name, rsmi):
                        is_hydrolysis = True
                        if reaction_name not in findings_json["atomic_checks"]["named_reactions"]:
                            findings_json["atomic_checks"]["named_reactions"].append(reaction_name)
                        break

                is_esterification = False
                for reaction_name in ESTERIFICATION_REACTIONS_OF_INTEREST:
                    if checker.check_reaction(reaction_name, rsmi):
                        is_esterification = True
                        if reaction_name not in findings_json["atomic_checks"]["named_reactions"]:
                            findings_json["atomic_checks"]["named_reactions"].append(reaction_name)
                        break

                if is_hydrolysis:
                    has_ester = False
                    for r in reactants:
                        if checker.check_fg("Ester", r):
                            has_ester = True
                            if "Ester" not in findings_json["atomic_checks"]["functional_groups"]:
                                findings_json["atomic_checks"]["functional_groups"].append("Ester")
                            break

                    has_acid = checker.check_fg("Carboxylic acid", product)
                    if has_acid:
                        if "Carboxylic acid" not in findings_json["atomic_checks"]["functional_groups"]:
                            findings_json["atomic_checks"]["functional_groups"].append("Carboxylic acid")

                    if has_ester and has_acid:
                        esterification_detected = True
                        # Find and add the specific co-occurrence constraint
                        for r_name in HYDROLYSIS_REACTIONS_OF_INTEREST:
                            if checker.check_reaction(r_name, rsmi):
                                constraint_targets = [r_name, "Ester", "Carboxylic acid"]
                                found_constraint = next((c for c in full_strategy_json["structural_constraints"] if c.get("type") == "co-occurrence" and c["details"].get("targets") == constraint_targets), None)
                                if found_constraint and found_constraint not in findings_json["structural_constraints"]:
                                    findings_json["structural_constraints"].append(copy.deepcopy(found_constraint))
                                break

                elif is_esterification:
                    has_acid = False
                    for r in reactants:
                        if checker.check_fg("Carboxylic acid", r):
                            has_acid = True
                            if "Carboxylic acid" not in findings_json["atomic_checks"]["functional_groups"]:
                                findings_json["atomic_checks"]["functional_groups"].append("Carboxylic acid")
                            break

                    has_ester = checker.check_fg("Ester", product)
                    if has_ester:
                        if "Ester" not in findings_json["atomic_checks"]["functional_groups"]:
                            findings_json["atomic_checks"]["functional_groups"].append("Ester")

                    if has_acid and has_ester:
                        esterification_detected = True
                        # Find and add the specific co-occurrence constraint
                        for r_name in ESTERIFICATION_REACTIONS_OF_INTEREST:
                            if checker.check_reaction(r_name, rsmi):
                                constraint_targets = [r_name, "Carboxylic acid", "Ester"]
                                found_constraint = next((c for c in full_strategy_json["structural_constraints"] if c.get("type") == "co-occurrence" and c["details"].get("targets") == constraint_targets), None)
                                if found_constraint and found_constraint not in findings_json["structural_constraints"]:
                                    findings_json["structural_constraints"].append(copy.deepcopy(found_constraint))
                                break

            except Exception:
                pass

        # Process children
        for child in node.get("children", []):
            new_depth = depth
            if node["type"] != "reaction": # If current node is chemical, depth increases
                new_depth = depth + 1
            dfs_traverse(child, new_depth)

    # Start traversal
    dfs_traverse(route)
    return esterification_detected, findings_json
