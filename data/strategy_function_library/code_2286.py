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


AMIDE_REDUCTION_REACTIONS = [
    "Reduction of primary amides to amines",
    "Reduction of secondary amides to amines",
    "Reduction of tertiary amides to amines",
]

def main(route) -> Tuple[bool, Dict]:
    """
    Detects if the synthesis uses a late-stage amide reduction strategy.
    Specifically, it checks if the final reaction step (depth=1) is one of the following types: 'Reduction of primary amides to amines', 'Reduction of secondary amides to amines', or 'Reduction of tertiary amides to amines'.
    """
    print("Starting late_stage_amide_reduction_strategy analysis")

    findings_template = {
        "atomic_checks": {
            "named_reactions": [],
            "ring_systems": [],
            "functional_groups": []
        },
        "structural_constraints": []
    }
    findings_json = copy.deepcopy(findings_template)

    # Validate route structure
    if not isinstance(route, dict) or "type" not in route:
        print("Invalid route structure")
        return False, findings_json

    # Ensure the root node is a molecule (the target compound)
    if route["type"] != "mol":
        print(f"Root node is not a molecule, but {route['type']}")
        return False, findings_json

    amide_reduction_found = False

    def dfs_traverse(node, depth=0):
        nonlocal amide_reduction_found, findings_json

        print(f"Traversing node of type {node['type']} at depth {depth}")

        if node["type"] == "mol":
            print(f"Molecule node: {node.get('smiles', 'No SMILES')}")

        # Check if this is a reaction node
        if node["type"] == "reaction" and "metadata" in node and "mapped_reaction_smiles" in node["metadata"]:
            rsmi = node["metadata"]["mapped_reaction_smiles"]
            print(f"Reaction node at depth {depth}: {rsmi}")

            # Check if this is the final step (depth 1)
            if depth == 1:  # The reaction leading to the target molecule
                print(f"Checking final step reaction: {rsmi}")

                # Check if this is an amide reduction reaction using the checker functions
                for rxn_name in AMIDE_REDUCTION_REACTIONS:
                    if checker.check_reaction(rxn_name, rsmi):
                        print(f"Found amide reduction reaction in final step: {rsmi}")
                        amide_reduction_found = True
                        findings_json["atomic_checks"]["named_reactions"].append(rxn_name)
                        # Add the structural constraint if the condition is met
                        findings_json["structural_constraints"].append({
                            "type": "positional",
                            "details": {
                                "target": [
                                    "Reduction of primary amides to amines",
                                    "Reduction of secondary amides to amines",
                                    "Reduction of tertiary amides to amines"
                                ],
                                "position": "last_stage"
                            }
                        })
                        break # Found one, no need to check others

        # Determine the next depth based on the current node's type
        next_depth = depth
        if node["type"] != "reaction": # If current node is 'mol', depth increases when going to 'reaction'
            next_depth = depth + 1

        # Traverse children with updated depth
        for child in node.get("children", []):
            dfs_traverse(child, next_depth)

    # Start traversal from the root node (target molecule)
    dfs_traverse(route)
    print(f"Amide reduction found: {amide_reduction_found}")
    return amide_reduction_found, findings_json
