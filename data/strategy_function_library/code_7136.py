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
    "Schotten-Baumann to ester",
    "Transesterification",
    "O-alkylation of carboxylic acids with diazo compounds",
    "Oxidative esterification of primary alcohols",
]

def main(route) -> Tuple[bool, Dict]:
    """
    Checks for the formation of an ester in the late stages (final two steps) of a synthesis. It identifies specific named reactions, including 'Esterification of Carboxylic Acids', 'Schotten-Baumann to ester', 'Transesterification', 'O-alkylation of carboxylic acids with diazo compounds', and 'Oxidative esterification of primary alcohols'.
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

    late_stage_esterification = False

    def dfs_traverse(node, depth=0):
        nonlocal late_stage_esterification, findings_json

        if node["type"] == "reaction" and depth <= 2:  # Check reactions up to depth 2
            if "mapped_reaction_smiles" in node.get("metadata", {}):
                rsmi = node["metadata"]["mapped_reaction_smiles"]

                for reaction_name in ESTERIFICATION_REACTIONS:
                    if checker.check_reaction(reaction_name, rsmi):
                        late_stage_esterification = True
                        if reaction_name not in findings_json["atomic_checks"]["named_reactions"]:
                            findings_json["atomic_checks"]["named_reactions"].append(reaction_name)
                        # Add the structural constraint if it's met
                        # This constraint is met if any of the specified reactions are found in the late stage
                        if {"type": "positional", "details": {"targets": ["Esterification of Carboxylic Acids", "Schotten-Baumann to ester", "Transesterification", "O-alkylation of carboxylic acids with diazo compounds", "Oxidative esterification of primary alcohols"], "position": "late_stage"}} not in findings_json["structural_constraints"]:
                            findings_json["structural_constraints"].append({"type": "positional", "details": {"targets": ["Esterification of Carboxylic Acids", "Schotten-Baumann to ester", "Transesterification", "O-alkylation of carboxylic acids with diazo compounds", "Oxidative esterification of primary alcohols"], "position": "late_stage"}})
                        break # Reaction found, no need to check other names

        for child in node.get("children", []):
            # Continue traversal to check other branches, unless already found
            if late_stage_esterification:
                return
            
            # New depth calculation logic
            new_depth = depth
            if node["type"] != "reaction": # If current node is chemical, depth increases
                new_depth = depth + 1
            # If current node is reaction, depth remains the same
            
            dfs_traverse(child, new_depth)

    dfs_traverse(route)
    return late_stage_esterification, findings_json