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


# Refactoring for Enumeration: Isolate the list of reaction types
AMIDE_FORMATION_REACTIONS = [
    "Acylation of Nitrogen Nucleophiles by Carboxylic Acids",
    "Carboxylic acid with primary amine to amide",
    "Acyl chloride with primary amine to amide (Schotten-Baumann)",
    "Acyl chloride with secondary amine to amide",
    "Ester with primary amine to amide",
    "Ester with secondary amine to amide",
    "Acylation of Nitrogen Nucleophiles by Acyl/Thioacyl/Carbamoyl Halides and Analogs_N",
]

def main(route) -> Tuple[bool, Dict]:
    """
    This function detects if the synthesis involves a late-stage amide formation (in the final step). It specifically checks if the final reaction matches one of the following predefined types: 'Acylation of Nitrogen Nucleophiles by Carboxylic Acids', 'Carboxylic acid with primary amine to amide', 'Acyl chloride with primary amine to amide (Schotten-Baumann)', 'Acyl chloride with secondary amine to amide', 'Ester with primary amine to amide', 'Ester with secondary amine to amide', or 'Acylation of Nitrogen Nucleophiles by Acyl/Thioacyl/Carbamoyl Halides and Analogs_N'.
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

    # Initialize flags
    amide_formation_detected = False

    def is_final_step(node, parent_route, depth):
        """Check if this is the final step in the synthesis (first reaction in retrosynthesis)"""
        return (
            depth == 1
            and parent_route["type"] == "mol"
            and node in parent_route.get("children", [])
        )

    def dfs_traverse(node, parent=None, depth=0):
        nonlocal amide_formation_detected, findings_json

        # Process reaction nodes
        if node["type"] == "reaction":
            # Check if this is the final step (first reaction in retrosynthesis)
            if parent and is_final_step(node, parent, depth):
                # Get reaction SMILES
                if "metadata" in node and "mapped_reaction_smiles" in node["metadata"]:
                    rsmi = node["metadata"]["mapped_reaction_smiles"]

                    # Check for amide formation reactions
                    is_amide_formation = False
                    for reaction_type in AMIDE_FORMATION_REACTIONS:
                        if checker.check_reaction(reaction_type, rsmi):
                            is_amide_formation = True
                            findings_json["atomic_checks"]["named_reactions"].append(reaction_type)
                            break

                    if is_amide_formation:
                        amide_formation_detected = True
                        # Add the structural constraint when detected
                        findings_json["structural_constraints"].append({
                            "type": "positional",
                            "details": {
                                "target": AMIDE_FORMATION_REACTIONS,
                                "position": "last_stage"
                            }
                        })
                        return

        # Determine the new depth for the recursive call
        new_depth = depth
        if node["type"] != "reaction":  # If current node is chemical, depth increases
            new_depth = depth + 1

        # Traverse children
        for child in node.get("children", []):
            if amide_formation_detected:
                break
            dfs_traverse(child, node, new_depth)

    # Start traversal
    dfs_traverse(route)

    return amide_formation_detected, findings_json
