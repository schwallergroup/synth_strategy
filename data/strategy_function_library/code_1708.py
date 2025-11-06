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


BOC_DEPROTECTION_REACTIONS = [
    "Boc amine deprotection",
    "Boc amine deprotection of guanidine",
    "Boc amine deprotection to NH-NH2",
    "Tert-butyl deprotection of amine",
]

def main(route) -> Tuple[bool, Dict]:
    """
    Detects if the synthetic route involves a specific Boc deprotection reaction as the final step, as defined in the BOC_DEPROTECTION_REACTIONS list.
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

    print("Starting Boc deprotection strategy analysis...")

    if not route or "type" not in route:
        print("Invalid route structure")
        return False, findings_json

    final_deprotection = False

    def dfs_traverse(node, depth=0):
        nonlocal final_deprotection, findings_json

        print(f"Analyzing node at depth {depth}, type: {node['type']}")

        # Check if this is a reaction node at the first step (depth 1)
        if node["type"] == "reaction" and depth == 1 and len(node.get("children", [])) > 0:
            print(f"Examining final reaction node at depth {depth}")

            try:
                if "metadata" in node and "mapped_reaction_smiles" in node["metadata"]:
                    rsmi = node["metadata"]["mapped_reaction_smiles"]
                    print(f"Reaction SMILES: {rsmi}")

                    # Check if this is a Boc deprotection reaction using various reaction types
                    for rxn_type in BOC_DEPROTECTION_REACTIONS:
                        if checker.check_reaction(rxn_type, rsmi):
                            print(f"Detected {rxn_type} in final step: {rsmi}")
                            final_deprotection = True
                            findings_json["atomic_checks"]["named_reactions"].append(rxn_type)
                            # Add the structural constraint if a Boc deprotection reaction is found at the final step
                            findings_json["structural_constraints"].append({
                                "type": "positional",
                                "details": {
                                    "targets": [
                                        "Boc amine deprotection",
                                        "Boc amine deprotection of guanidine",
                                        "Boc amine deprotection to NH-NH2",
                                        "Tert-butyl deprotection of amine"
                                    ],
                                    "position": "last_stage"
                                }
                            })
                            return
                else:
                    print("No reaction SMILES found in metadata")
            except Exception as e:
                print(f"Error analyzing reaction at depth {depth}: {e}")

        # Traverse children
        for child in node.get("children", []):
            # New logic for depth calculation
            new_depth = depth
            if node['type'] != 'reaction': # If current node is chemical, depth increases
                new_depth = depth + 1
            # If current node is reaction, depth remains the same
            dfs_traverse(child, new_depth)

    dfs_traverse(route)
    print(f"Final result: {final_deprotection}")
    return final_deprotection, findings_json
