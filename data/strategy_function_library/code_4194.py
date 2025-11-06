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


GROUPS_TO_PRESERVE = ["Trifluoro group", "Nitrile"]

def main(route) -> Tuple[bool, Dict]:
    """
    Detects if a synthetic route maintains a specific set of functional groups, as defined in the `GROUPS_TO_PRESERVE` list (e.g., 'Trifluoro group', 'Nitrile'), throughout all reaction steps. The strategy is only considered valid if the final product initially contains all specified groups.
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

    result = True

    # First, check if the final product has all the required groups
    if route["type"] == "mol":
        all_groups_present_in_final_product = True
        for group in GROUPS_TO_PRESERVE:
            if checker.check_fg(group, route["smiles"]):
                findings_json["atomic_checks"]["functional_groups"].append(group)
            else:
                all_groups_present_in_final_product = False
                result = False # If final product is missing any group, no need to check further
        if all_groups_present_in_final_product:
            findings_json["structural_constraints"].append({
                "type": "co-occurrence",
                "details": {
                    "targets": [
                        "Trifluoro group",
                        "Nitrile"
                    ],
                    "scope": "final_product"
                }
            })
    else:
        result = False  # Root must be a molecule

    if not result: # If initial check failed, no need to proceed with traversal
        return result, findings_json

    # Track if groups are preserved when present
    preservation_status = {group: True for group in GROUPS_TO_PRESERVE}

    def dfs_traverse(node, depth=0):
        nonlocal preservation_status, result, findings_json

        # Early exit if any group has already been found to be lost
        if not all(preservation_status.values()):
            return

        if node["type"] == "reaction" and "metadata" in node and "mapped_reaction_smiles" in node["metadata"]:
            try:
                rsmi = node["metadata"]["mapped_reaction_smiles"]
                reactants = rsmi.split(">")[0].split(".")
                product = rsmi.split(">")[-1]

                for group in GROUPS_TO_PRESERVE:
                    # Check if a group present in any reactant is lost in the product
                    group_present_in_reactants = any(checker.check_fg(group, r) for r in reactants)
                    group_present_in_product = checker.check_fg(group, product)

                    if group_present_in_reactants and not group_present_in_product:
                        preservation_status[group] = False
                        result = False # A group was lost, so the overall result is False
                        # No specific atomic check for 'loss', but the structural constraint will capture it.

            except Exception:
                # Silently ignore reactions that fail to parse
                pass

        # Traverse children
        for child in node.get("children", []):
            # New logic for depth calculation
            if node["type"] == "reaction":
                # If current node is a reaction, depth remains the same for children
                dfs_traverse(child, depth)
            else:
                # If current node is not a reaction (e.g., chemical), depth increases
                dfs_traverse(child, depth + 1)

    # Start traversal from the root
    dfs_traverse(route)

    # After traversal, add structural constraints based on preservation status
    for group, preserved in preservation_status.items():
        if not preserved:
            # If a group was NOT preserved, it means the negation constraint was violated.
            # We record the violation by NOT adding the negation constraint to findings_json.
            # However, the prompt asks to add 'detected' elements. If the overall result is False
            # due to loss, we don't add the 'negation' constraint as it wasn't 'met'.
            # The 'result' variable already captures this.
            pass # The 'result' being False is the indicator.
        else:
            # If a group was preserved, add its corresponding negation constraint as 'met'
            findings_json["structural_constraints"].append({
                "type": "negation",
                "details": {
                    "target": f"loss_of_{group.replace(' ', '_')}",
                    "scope": "all_steps"
                }
            })

    # Return true if all groups were preserved throughout the synthesis
    return result, findings_json
