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

def main(route) -> Tuple[bool, Dict]:
    """
    Detects the use of specific, known named reactions for the protection of a carboxylic acid as an ester or the deprotection of an ester to a carboxylic acid.
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

    has_acid_protection_or_deprotection = False

    # Define the structural constraint object from the original JSON for easy access
    structural_constraint_obj = {
      "type": "count",
      "details": {
        "targets": [
          "Protection of carboxylic acid",
          "Esterification of Carboxylic Acids",
          "O-alkylation of carboxylic acids with diazo compounds",
          "Ester saponification (methyl deprotection)",
          "Ester saponification (alkyl deprotection)",
          "COOH ethyl deprotection",
          "Deprotection of carboxylic acid",
          "Hydrolysis or Hydrogenolysis of Carboxylic Esters or Thioesters"
        ],
        "operator": ">=",
        "value": 1
      }
    }

    def dfs_traverse(node, depth=0):
        nonlocal has_acid_protection_or_deprotection, findings_json

        if node["type"] == "reaction" and "metadata" in node and "mapped_reaction_smiles" in node["metadata"]:
            rsmi = node["metadata"]["mapped_reaction_smiles"]

            # Check if this is a known carboxylic acid protection reaction
            protection_reactions = [
                "Protection of carboxylic acid",
                "Esterification of Carboxylic Acids",
                "O-alkylation of carboxylic acids with diazo compounds"
            ]
            for r_name in protection_reactions:
                if checker.check_reaction(r_name, rsmi):
                    has_acid_protection_or_deprotection = True
                    if r_name not in findings_json["atomic_checks"]["named_reactions"]:
                        findings_json["atomic_checks"]["named_reactions"].append(r_name)
                    # No return here, continue to check other reactions in case multiple apply

            # Check if this is a known ester deprotection reaction
            deprotection_reactions = [
                "Ester saponification (methyl deprotection)",
                "Ester saponification (alkyl deprotection)",
                "COOH ethyl deprotection",
                "Deprotection of carboxylic acid",
                "Hydrolysis or Hydrogenolysis of Carboxylic Esters or Thioesters"
            ]
            for r_name in deprotection_reactions:
                if checker.check_reaction(r_name, rsmi):
                    has_acid_protection_or_deprotection = True
                    if r_name not in findings_json["atomic_checks"]["named_reactions"]:
                        findings_json["atomic_checks"]["named_reactions"].append(r_name)
                    # No return here, continue to check other reactions in case multiple apply

        # Traverse children
        for child in node.get("children", []):
            new_depth = depth
            if node["type"] != "reaction": # If current node is not a reaction (e.g., chemical), increase depth
                new_depth = depth + 1
            dfs_traverse(child, new_depth)

    dfs_traverse(route)

    # Add the structural constraint if any relevant reaction was found
    if has_acid_protection_or_deprotection and structural_constraint_obj not in findings_json["structural_constraints"]:
        findings_json["structural_constraints"].append(structural_constraint_obj)

    return has_acid_protection_or_deprotection, findings_json
