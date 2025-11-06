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


LATE_STAGE_AMIDE_FORMATION_REACTIONS = [
    "Acylation of Nitrogen Nucleophiles by Acyl/Thioacyl/Carbamoyl Halides and Analogs_N",
    "Acylation of Nitrogen Nucleophiles by Acyl/Thioacyl/Carbamoyl Halides and Analogs_OS",
    "Acylation of Nitrogen Nucleophiles by Carboxylic Acids",
    "Acyl chloride with ammonia to amide",
    "Acyl chloride with primary amine to amide (Schotten-Baumann)",
    "Acyl chloride with secondary amine to amide",
    "Carboxylic acid with primary amine to amide",
    "Ester with ammonia to amide",
    "Ester with primary amine to amide",
    "Ester with secondary amine to amide",
    "Schotten-Baumann_amide",
]

def main(route) -> Tuple[bool, Dict]:
    """
    Checks for the formation of an amide bond in the final synthetic step, identified from a specific list of named reactions including various acylations (e.g., with acyl halides, carboxylic acids, esters) and Schotten-Baumann type reactions.
    """
    final_step_is_amide_formation = False

    findings_template = {
        "atomic_checks": {
            "named_reactions": [],
            "ring_systems": [],
            "functional_groups": []
        },
        "structural_constraints": []
    }
    findings_json = copy.deepcopy(findings_template)

    def dfs_traverse(node, depth=0):
        nonlocal final_step_is_amide_formation, findings_json

        # Check if this is the target molecule (depth 0) and it has a reaction child
        if node["type"] == "mol" and depth == 0:
            if node.get("children", []) and node["children"][0]["type"] == "reaction":
                reaction_node = node["children"][0]
                if "metadata" in reaction_node and "mapped_reaction_smiles" in reaction_node["metadata"]:
                    rsmi = reaction_node["metadata"]["mapped_reaction_smiles"]

                    # Check for amide formation reactions using the checker function
                    for reaction_type in LATE_STAGE_AMIDE_FORMATION_REACTIONS:
                        if checker.check_reaction(reaction_type, rsmi):
                            final_step_is_amide_formation = True
                            findings_json["atomic_checks"]["named_reactions"].append(reaction_type)
                            # Add structural constraint for last_stage reaction
                            findings_json["structural_constraints"].append({
                                "type": "positional",
                                "details": {
                                    "target": LATE_STAGE_AMIDE_FORMATION_REACTIONS,
                                    "position": "last_stage"
                                }
                            })
                            return

        for child in node.get("children", []):
            # Determine the new depth based on the current node's type
            new_depth = depth
            if node["type"] != "reaction": # If current node is chemical, depth increases
                new_depth = depth + 1
            # If current node is reaction, depth remains the same
            
            dfs_traverse(child, new_depth)

    # Ensure we start with a molecule node
    if route["type"] == "mol":
        dfs_traverse(route)
    else:
        print("Root node is not a molecule")

    return final_step_is_amide_formation, findings_json
