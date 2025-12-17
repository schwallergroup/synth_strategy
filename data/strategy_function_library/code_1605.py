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


REDUCTIVE_AMINATION_REACTIONS = [
    "Reductive amination with aldehyde",
    "Reductive amination with ketone",
    "Reductive amination with alcohol",
    "Mignonac reaction",
]

def main(route) -> Tuple[bool, Dict]:
    """
    Detects late-stage C-N bond formation by checking for a specific set of named reactions: 'Reductive amination with aldehyde', 'Reductive amination with ketone', 'Reductive amination with alcohol', and 'Mignonac reaction'. The search is limited to the final three steps of the synthesis.
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

    reductive_amination_detected = False

    def dfs_traverse(node, depth=0):
        nonlocal reductive_amination_detected, findings_json

        if (
            node["type"] == "reaction" and depth <= 2
        ):  # Check late-stage reactions
            if "metadata" in node and "mapped_reaction_smiles" in node["metadata"]:
                rsmi = node["metadata"]["mapped_reaction_smiles"]

                # Check if this is a reductive amination reaction using a predefined list
                is_reductive_amination = False
                for name in REDUCTIVE_AMINATION_REACTIONS:
                    if checker.check_reaction(name, rsmi):
                        is_reductive_amination = True
                        findings_json["atomic_checks"]["named_reactions"].append(name)
                        break

                if is_reductive_amination:
                    reductive_amination_detected = True
                    # Add structural constraint if detected within the last three stages
                    findings_json["structural_constraints"].append({
                        "type": "positional",
                        "details": {
                            "target": [
                                "Reductive amination with aldehyde",
                                "Reductive amination with ketone",
                                "Reductive amination with alcohol",
                                "Mignonac reaction"
                            ],
                            "position": "last_three_stages"
                        }
                    })

        for child in node.get("children", []):
            # New logic for depth calculation
            new_depth = depth
            if node["type"] != "reaction": # If current node is chemical, depth increases
                new_depth = depth + 1
            # If current node is reaction, depth remains the same
            
            dfs_traverse(child, new_depth)

    dfs_traverse(route)
    return reductive_amination_detected, findings_json
