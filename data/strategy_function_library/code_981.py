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

def main(route) -> Tuple[bool, Dict]:
    """
    This function detects if the synthetic route involves a sequence of
    nitro group introduction followed by reduction to amine.
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

    nitration_steps = []
    reduction_steps = []
    result = False

    def dfs_traverse(node, depth=0):
        nonlocal nitration_steps, reduction_steps, findings_json
        if node["type"] == "reaction":
            try:
                rsmi = node["metadata"]["mapped_reaction_smiles"]

                # Check for nitration reactions
                nitration_reactions = [
                    "Aromatic nitration with HNO3",
                    "Aromatic nitration with NO3 salt",
                    "Aromatic nitration with NO2 salt",
                    "Aromatic nitration with alkyl NO2",
                    "Non-aromatic nitration with HNO3"
                ]
                for r_name in nitration_reactions:
                    if checker.check_reaction(r_name, rsmi):
                        nitration_steps.append(depth)
                        if r_name not in findings_json["atomic_checks"]["named_reactions"]:
                            findings_json["atomic_checks"]["named_reactions"].append(r_name)
                        print(f"Detected nitration at depth {depth}")
                        break # Only add one nitration type per reaction node

                # Check for nitro reduction to amine
                reduction_reaction_name = "Reduction of nitro groups to amines"
                if checker.check_reaction(reduction_reaction_name, rsmi):
                    reduction_steps.append(depth)
                    if reduction_reaction_name not in findings_json["atomic_checks"]["named_reactions"]:
                        findings_json["atomic_checks"]["named_reactions"].append(reduction_reaction_name)
                    print(f"Detected nitro reduction at depth {depth}")
            except Exception as e:
                print(f"Error processing reaction: {e}")

        # Determine the next depth based on the current node's type
        next_depth = depth
        if node["type"] != "reaction": # If current node is chemical, depth increases
            next_depth = depth + 1

        # Traverse children
        for child in node.get("children", []):
            dfs_traverse(child, next_depth)

    # Start traversal
    dfs_traverse(route)

    # Check if there's a nitration followed by reduction
    for nitration_depth in nitration_steps:
        for reduction_depth in reduction_steps:
            # In retrosynthetic traversal, lower depth = later in synthesis (earlier in retrosynthesis)
            if (
                reduction_depth < nitration_depth
            ):  # Check if reduction happens after nitration in forward synthesis
                print("Detected nitro-reduction-amination sequence")
                result = True
                # Add the structural constraint to findings_json
                structural_constraint_obj = {
                    "type": "sequence",
                    "details": {
                        "before_event": [
                            "Aromatic nitration with HNO3",
                            "Aromatic nitration with NO3 salt",
                            "Aromatic nitration with NO2 salt",
                            "Aromatic nitration with alkyl NO2",
                            "Non-aromatic nitration with HNO3"
                        ],
                        "after_event": [
                            "Reduction of nitro groups to amines"
                        ]
                    }
                }
                if structural_constraint_obj not in findings_json["structural_constraints"]:
                    findings_json["structural_constraints"].append(structural_constraint_obj)
                # No need to break here, as we want to find all such sequences if multiple exist
                # The problem statement implies returning as soon as one is found, but for findings, we collect all.
                # However, the original function returns True and exits on first match, so we'll keep that behavior for 'result'.
                return result, findings_json # Exit early as per original function's behavior

    return result, findings_json