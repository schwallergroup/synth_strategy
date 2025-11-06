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
    Detects late-stage (depth <= 1) ester hydrolysis by checking if a reaction matches any of the reaction types defined in the `ESTER_HYDROLYSIS_REACTIONS` list, which includes various saponification and general hydrolysis reactions.
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

    ester_hydrolysis_found = False
    min_depth_with_ester_hydrolysis = float("inf")

    def dfs_traverse(node, current_depth=0):
        nonlocal ester_hydrolysis_found, min_depth_with_ester_hydrolysis, findings_json

        if node["type"] == "reaction":
            # Get depth from metadata if available, otherwise use current_depth
            depth = node.get("metadata", {}).get("depth", current_depth)

            # Extract reactants and product
            try:
                rsmi = node["metadata"].get("rsmi", "")
                if not rsmi:
                    return

                # Check if this is an ester hydrolysis reaction using multiple reaction types
                is_ester_hydrolysis = False
                for rxn_name in ESTER_HYDROLYSIS_REACTIONS:
                    if checker.check_reaction(rxn_name, rsmi):
                        is_ester_hydrolysis = True
                        findings_json["atomic_checks"]["named_reactions"].append(rxn_name)
                        break

                if is_ester_hydrolysis:
                    ester_hydrolysis_found = True

                    # Keep track of the minimum depth where we found ester hydrolysis
                    if depth < min_depth_with_ester_hydrolysis:
                        min_depth_with_ester_hydrolysis = depth
            except Exception as e:
                print(f"Error processing reaction node: {e}")

        # Determine the next depth based on the current node's type
        next_depth = current_depth
        if node["type"] != "reaction": # If current node is chemical, depth increases
            next_depth = current_depth + 1

        # Traverse children with the determined next depth
        for child in node.get("children", []):
            dfs_traverse(child, next_depth)

    # Start traversal from root
    dfs_traverse(route)

    # Consider it late-stage if it's at depth 0 or 1
    is_late_stage = min_depth_with_ester_hydrolysis <= 1

    result = is_late_stage and ester_hydrolysis_found

    if result:
        # Add the structural constraint if the overall condition is met
        findings_json["structural_constraints"].append({
            "type": "positional",
            "details": {
                "target": "ester_hydrolysis",
                "position": "late_stage (depth <= 1)"
            }
        })

    return result, findings_json
