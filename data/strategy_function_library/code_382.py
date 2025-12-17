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


# This list should be defined at the module level, outside the main function.
AMIDE_FORMATION_REACTIONS = [
    "Acylation of Nitrogen Nucleophiles by Acyl/Thioacyl/Carbamoyl Halides and Analogs_N",
    "Acylation of Nitrogen Nucleophiles by Carboxylic Acids",
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
    This function detects a synthetic strategy involving late-stage amide formation
    as the final step in a linear synthesis.
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

    amide_formation_at_final_step = False
    is_linear_synthesis = True
    reaction_count = 0

    def dfs_traverse(node, depth=0, branch_id=0):
        nonlocal amide_formation_at_final_step, is_linear_synthesis, reaction_count, findings_json

        # Check if this is a reaction node
        if node["type"] == "reaction" and "metadata" in node and "mapped_reaction_smiles" in node["metadata"]:
            reaction_count += 1
            rsmi = node["metadata"]["mapped_reaction_smiles"]

            # In retrosynthesis, depth 1 is the final step in forward synthesis
            is_final_step = depth == 1

            if is_final_step:
                # Check for amide formation reactions using the checker function
                is_amide_formation = False
                for rxn_type in AMIDE_FORMATION_REACTIONS:
                    if checker.check_reaction(rxn_type, rsmi):
                        is_amide_formation = True
                        findings_json["atomic_checks"]["named_reactions"].append(rxn_type)
                        break

                if is_amide_formation:
                    amide_formation_at_final_step = True

        # Check for branching in the synthesis route
        if node["type"] == "mol" and len(node.get("children", [])) > 1:
            # More than one reaction from this molecule indicates branching
            is_linear_synthesis = False

        # Determine the depth for the recursive call based on the current node's type
        next_depth = depth
        if node["type"] != "reaction": # This means it's a 'mol' node, or any other type that should increase depth
            next_depth = depth + 1

        # Continue traversing
        child_count = 0
        for child in node.get("children", []):
            child_count += 1
            # Pass a new branch_id for each child of a molecule with multiple children
            new_branch_id = branch_id
            if node["type"] == "mol" and child_count > 1:
                new_branch_id = branch_id + 1
            dfs_traverse(child, next_depth, new_branch_id)

    # Start traversal from the root
    dfs_traverse(route)

    result = amide_formation_at_final_step and is_linear_synthesis

    if amide_formation_at_final_step:
        findings_json["structural_constraints"].append({
            "type": "positional",
            "details": {
                "target": "amide_formation",
                "position": "last_stage"
            }
        })
    if is_linear_synthesis:
        # This constraint is met if is_linear_synthesis remains True
        # The original JSON defines this as a negation, so if it's true, the negation is met.
        findings_json["structural_constraints"].append({
            "type": "negation",
            "details": {
                "target": "route_branching"
            }
        })

    return result, findings_json
