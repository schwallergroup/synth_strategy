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


BOC_DEPROTECTION_REACTION_TYPES = [
    "Boc amine deprotection",
    "Boc amine deprotection of guanidine",
    "Boc amine deprotection to NH-NH2",
    "Tert-butyl deprotection of amine",
]

def main(route) -> Tuple[bool, Dict]:
    """
    Detects if a final synthetic step is a specific Boc-related deprotection, as defined in the BOC_DEPROTECTION_REACTION_TYPES list. The check is performed on all leaf reaction nodes in the synthesis tree.
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

    # Helper function to check if a reaction is a Boc deprotection
    def is_boc_deprotection(reaction_node):
        nonlocal findings_json
        if not reaction_node.get("metadata", {}).get("rsmi"):
            return False

        rsmi = reaction_node["metadata"]["mapped_reaction_smiles"]

        # Check if the reaction is a Boc deprotection using the checker function
        for boc_type in BOC_DEPROTECTION_REACTION_TYPES:
            if checker.check_reaction(boc_type, rsmi):
                print(f"Found {boc_type} as final step: {rsmi}")
                findings_json["atomic_checks"]["named_reactions"].append(boc_type)
                # Add the corresponding structural constraint
                findings_json["structural_constraints"].append({
                    "type": "positional",
                    "details": {
                        "target": boc_type,
                        "position": "last_stage"
                    }
                })
                return True

        return False

    # Find all leaf reaction nodes (final steps in forward synthesis)
    found_boc_deprotection = [False]  # Using list to allow modification in nested function

    def find_leaf_reactions(node, is_child_of_reaction=False):
        nonlocal found_boc_deprotection
        if found_boc_deprotection[0]: # Early exit
            return

        if node["type"] == "reaction":
            # Check if this reaction node has only molecule children (leaf reaction)
            all_children_are_mols = all(
                child["type"] == "mol" for child in node.get("children", [])
            )

            if all_children_are_mols:
                # This is a leaf reaction node (final step in forward synthesis)
                if is_boc_deprotection(node):
                    found_boc_deprotection[0] = True
                    return

            # Continue traversal for non-leaf reaction nodes
            for child in node.get("children", []):
                find_leaf_reactions(child, True)

        elif node["type"] == "mol" and not is_child_of_reaction:
            # For molecule nodes that aren't children of reaction nodes
            for child in node.get("children", []):
                find_leaf_reactions(child, False)

    # Start traversal from the root
    find_leaf_reactions(route)

    return found_boc_deprotection[0], findings_json


def dfs_traverse(node, depth=0, visited=None):
    if visited is None:
        visited = set()

    node_id = id(node)
    if node_id in visited:
        return
    visited.add(node_id)

    node["depth"] = depth

    for child in node.get("children", []):
        new_depth = depth
        if node["type"] != "reaction": # This means current node is 'chemical' or 'mol'
            new_depth = depth + 1
        # If current node is 'reaction', new_depth remains 'depth'
        dfs_traverse(child, new_depth, visited)
