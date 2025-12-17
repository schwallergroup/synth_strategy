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


AROMATIC_HALOGENATION_REACTIONS = [
    "Aromatic fluorination",
    "Aromatic chlorination",
    "Aromatic bromination",
    "Aromatic iodination",
]

def main(route) -> Tuple[bool, Dict]:
    """
    Detects if the synthesis route involves a late-stage aromatic halogenation in the final step.
    This is determined by checking if the final reaction is one of the types specified in the
    AROMATIC_HALOGENATION_REACTIONS list.
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

    halogenation_detected = False

    def dfs_traverse(node, depth=0):
        nonlocal halogenation_detected, findings_json

        if node["type"] == "reaction" and node.get("metadata", {}).get("rsmi"):
            # The final reaction step is at depth=1.
            if depth == 1:
                rsmi = node["metadata"]["mapped_reaction_smiles"]

                is_halogenation_found_in_step = False
                for reaction_name in AROMATIC_HALOGENATION_REACTIONS:
                    if checker.check_reaction(reaction_name, rsmi):
                        is_halogenation_found_in_step = True
                        if reaction_name not in findings_json["atomic_checks"]["named_reactions"]:
                            findings_json["atomic_checks"]["named_reactions"].append(reaction_name)

                if is_halogenation_found_in_step:
                    halogenation_detected = True
                    # Add the structural constraint if the condition is met
                    # This corresponds to the 'positional' constraint in the strategy JSON
                    findings_json["structural_constraints"].append({
                        "type": "positional",
                        "details": {
                            "targets": [
                                "Aromatic fluorination",
                                "Aromatic chlorination",
                                "Aromatic bromination",
                                "Aromatic iodination"
                            ],
                            "position": "last_stage"
                        }
                    })

        # Continue traversal
        for child in node.get("children", []):
            # New logic for depth calculation:
            # Depth increases only when traversing from a chemical node to a reaction node.
            # Depth remains the same when traversing from a reaction node to a chemical node.
            if node["type"] == "reaction":
                dfs_traverse(child, depth)
            else: # node['type'] == 'chemical'
                dfs_traverse(child, depth + 1)

    # Start traversal from the root
    dfs_traverse(route)
    return halogenation_detected, findings_json
