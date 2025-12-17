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

def main(route) -> Tuple[bool, Dict]:
    """
    This function detects a synthetic strategy involving nitrile hydrolysis to carboxylic acid
    in the late stage of synthesis.
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

    nitrile_hydrolysis_found = False
    max_depth = 10  # Arbitrary large number
    late_stage_threshold = 3  # Consider reactions within depth 3 as late stage

    def dfs_traverse(node, depth=0):
        nonlocal nitrile_hydrolysis_found, findings_json

        if node["type"] == "reaction" and "metadata" in node and "mapped_reaction_smiles" in node["metadata"]:
            rsmi = node["metadata"]["mapped_reaction_smiles"]

            try:
                # Check if this is a nitrile to carboxylic acid oxidation reaction
                if checker.check_reaction("Oxidation of nitrile to carboxylic acid", rsmi):
                    findings_json["atomic_checks"]["named_reactions"].append("Oxidation of nitrile to carboxylic acid")
                    # Check if this is in the late stage (closer to final product)
                    if depth <= late_stage_threshold:
                        nitrile_hydrolysis_found = True
                        # Add the structural constraint if both conditions are met
                        findings_json["structural_constraints"].append({
                            "type": "positional",
                            "details": {
                                "target": "Oxidation of nitrile to carboxylic acid",
                                "position": "late_stage"
                            }
                        })
            except Exception as e:
                print(f"Error processing reaction: {e}")

        # Traverse children
        for child in node.get("children", []):
            new_depth = depth
            if node["type"] != "reaction": # If current node is chemical, depth increases
                new_depth = depth + 1
            dfs_traverse(child, new_depth)

    # Start traversal from the root
    dfs_traverse(route)

    print(f"Total synthesis depth: {max_depth}")
    print(f"Late stage threshold: {late_stage_threshold}")
    print(f"Nitrile hydrolysis in late stage: {nitrile_hydrolysis_found}")

    return nitrile_hydrolysis_found, findings_json
