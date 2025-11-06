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


REDUCTIVE_AMINATION_VARIANTS = [
    "Reductive amination with aldehyde",
    "Reductive amination with ketone",
    "Reductive amination with alcohol",
    "Mignonac reaction",
]

def main(route) -> Tuple[bool, Dict]:
    """
    Detects a late-stage (synthesis steps 1 or 2) C-N bond formation using one of a specific list of reductive amination variants: 'Reductive amination with aldehyde', 'Reductive amination with ketone', 'Reductive amination with alcohol', and 'Mignonac reaction'.
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

    has_reductive_amination = False

    def dfs_traverse(node, depth=0):
        nonlocal has_reductive_amination, findings_json

        # Check for reductive amination in late-stage reactions (depth ≤ 2)
        if node["type"] == "reaction" and depth <= 2:
            if "mapped_reaction_smiles" in node["metadata"]:
                rsmi = node["metadata"]["mapped_reaction_smiles"]

                # Use checker function to identify reductive amination reactions
                for reaction_name in REDUCTIVE_AMINATION_VARIANTS:
                    if checker.check_reaction(reaction_name, rsmi):
                        has_reductive_amination = True
                        if reaction_name not in findings_json["atomic_checks"]["named_reactions"]:
                            findings_json["atomic_checks"]["named_reactions"].append(reaction_name)
                        # Add structural constraint if a relevant reaction is found at the correct depth
                        # This assumes that finding any of these reactions at depth <= 2 satisfies the positional constraint
                        positional_constraint_obj = {
                            "type": "positional",
                            "details": {
                                "target": [
                                    "Reductive amination with aldehyde",
                                    "Reductive amination with ketone",
                                    "Reductive amination with alcohol",
                                    "Mignonac reaction"
                                ],
                                "position": "within_first_2_retrosynthetic_steps"
                            }
                        }
                        if positional_constraint_obj not in findings_json["structural_constraints"]:
                            findings_json["structural_constraints"].append(positional_constraint_obj)
                        break # Exit loop once a match is found for this reaction node

        # Continue traversing the tree
        for child in node.get("children", []):
            # New logic: depth increases only when going from chemical to reaction
            if node["type"] == "reaction":
                dfs_traverse(child, depth)
            else: # node["type"] == "chemical" or any other type
                dfs_traverse(child, depth + 1)

    # Start traversal from the root
    if route["type"] == "mol":
        dfs_traverse(route)
    else:
        dfs_traverse(route)

    return has_reductive_amination, findings_json
