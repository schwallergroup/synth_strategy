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


ESTER_HYDROLYSIS_REACTIONS = [
    "Hydrolysis or Hydrogenolysis of Carboxylic Esters or Thioesters",
    "Ester saponification (methyl deprotection)",
    "Ester saponification (alkyl deprotection)",
]

def main(route) -> Tuple[bool, Dict]:
    """
    Detects late-stage esterification by identifying its retrosynthetic equivalent—ester hydrolysis—in the final synthetic step.
    This is achieved by checking for specific reaction types defined in ESTER_HYDROLYSIS_REACTIONS, including general hydrolysis, hydrogenolysis, and saponification of esters.
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

    esterification_at_final_step = False

    def dfs_traverse(node, depth=0):
        nonlocal esterification_at_final_step, findings_json

        # Set depth in the node for reference
        node["depth"] = depth

        if node["type"] == "reaction" and "metadata" in node and "mapped_reaction_smiles" in node["metadata"]:
            rsmi = node["metadata"]["mapped_reaction_smiles"]

            # Check if this is a late-stage step (depth 0 or 1)
            if depth <= 1:
                for reaction_name in ESTER_HYDROLYSIS_REACTIONS:
                    if checker.check_reaction(reaction_name, rsmi):
                        esterification_at_final_step = True
                        findings_json["atomic_checks"]["named_reactions"].append(reaction_name)
                        # Add the structural constraint if the condition is met
                        findings_json["structural_constraints"].append({
                            "type": "positional",
                            "details": {
                                "targets": [
                                    "Hydrolysis or Hydrogenolysis of Carboxylic Esters or Thioesters",
                                    "Ester saponification (methyl deprotection)",
                                    "Ester saponification (alkyl deprotection)"
                                ],
                                "target_quantifier": "any",
                                "position": "late_stage (depth <= 1)"
                            }
                        })
                        break  # A match is found, no need to check other reaction types

        # Traverse children with incremented depth
        for child in node.get("children", []):
            # New logic: depth increases only when traversing from chemical to reaction
            if node["type"] == "reaction":
                dfs_traverse(child, depth)
            else: # node['type'] == 'chemical'
                dfs_traverse(child, depth + 1)

    # Start traversal from the root
    try:
        dfs_traverse(route)
    except Exception as e:
        print(f"Error during traversal: {e}")

    return esterification_at_final_step, findings_json
