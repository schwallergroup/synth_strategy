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


DECARBOXYLATION_REACTION_TYPES = [
    "Decarboxylation",
    "Ketonization by decarboxylation of carbonic acids",
    "Ketonization by decarboxylation of acid halides",
    "decarboxylative_coupling",
]

def main(route) -> Tuple[bool, Dict]:
    """
    Checks if the final synthetic step is a decarboxylation-type reaction. The check is performed by matching the reaction against a defined list of named decarboxylation reactions (see DECARBOXYLATION_REACTION_TYPES).
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

    # Track if we found decarboxylation at late stage
    found_late_decarboxylation = False

    def calculate_depth(node, current_node):
        """Calculate depth of a reaction node in the synthesis route"""
        if node == current_node:
            return 0

        for child in node.get("children", []):
            if child == current_node:
                return 1

            for grandchild in child.get("children", []):
                if grandchild == current_node:
                    return 2

        return 3  # Depth > 2 (not late-stage)

    def dfs_traverse(node, depth=0):
        nonlocal found_late_decarboxylation, findings_json

        if node["type"] == "reaction":
            # Calculate depth based on position in the tree
            # The original calculate_depth function is not used here for the new logic
            # We use the 'depth' parameter passed to dfs_traverse directly.

            # Only check the final reaction step (depth=1)
            if depth <= 1:
                try:
                    rsmi = node["metadata"]["mapped_reaction_smiles"]

                    # Check for decarboxylation reaction types
                    for rxn_type in DECARBOXYLATION_REACTION_TYPES:
                        if checker.check_reaction(rxn_type, rsmi):
                            found_late_decarboxylation = True
                            findings_json["atomic_checks"]["named_reactions"].append(rxn_type)
                            # If a late-stage decarboxylation is found, add the structural constraint
                            findings_json["structural_constraints"].append({
                                "type": "positional",
                                "details": {
                                    "targets": [
                                        "Decarboxylation",
                                        "Ketonization by decarboxylation of carbonic acids",
                                        "Ketonization by decarboxylation of acid halides",
                                        "decarboxylative_coupling"
                                    ],
                                    "position": "last_stage"
                                }
                            })
                            return
                except Exception as e:
                    # Error analyzing reaction, pass silently
                    pass

        # Traverse children
        for child in node.get("children", []):
            # New depth calculation logic:
            # Depth increases only when traversing from a chemical node to a reaction node.
            # Depth remains the same when traversing from a reaction node to a chemical node.
            if node["type"] == "reaction":
                # From reaction to chemical, depth remains the same
                dfs_traverse(child, depth)
            else:
                # From chemical to reaction, depth increases
                dfs_traverse(child, depth + 1)

    # Start traversal
    # The initial call to dfs_traverse should start with depth 0 for the root node.
    dfs_traverse(route, depth=0)

    return found_late_decarboxylation, findings_json
