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

def main(route) -> Tuple[bool, Dict]:
    """
    Detects if a trifluoromethyl group is introduced early in the synthesis.
    Early introduction means the CF3 group is added in the first half of the synthesis steps.
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

    cf3_introduction_depths = []
    max_depth = 0
    target_has_cf3 = False
    result = False

    # First check if the target molecule contains a CF3 group
    if route["type"] == "mol":
        target_has_cf3 = checker.check_fg("Trifluoro group", route["smiles"])
        if target_has_cf3:
            findings_json["atomic_checks"]["functional_groups"].append("Trifluoro group")

        if not target_has_cf3:
            return False, findings_json

    def dfs_traverse(node, depth=0):
        nonlocal cf3_introduction_depths, max_depth, findings_json

        max_depth = max(max_depth, depth)

        if node["type"] == "reaction" and "metadata" in node and "mapped_reaction_smiles" in node["metadata"]:
            try:
                rsmi = node["metadata"]["mapped_reaction_smiles"]
                reactants_smiles = rsmi.split(">")[0].split(".")
                product_smiles = rsmi.split(">")[-1]

                # Check if CF3 is in the product
                product_has_cf3 = checker.check_fg("Trifluoro group", product_smiles)

                # Check if CF3 is in any reactant
                reactants_have_cf3 = any(
                    checker.check_fg("Trifluoro group", reactant) for reactant in reactants_smiles
                )

                # If CF3 is in product but not in reactants, it was introduced in this reaction
                if product_has_cf3 and not reactants_have_cf3:
                    cf3_introduction_depths.append(depth)
                    # Record the introduction of Trifluoro group
                    if "Trifluoro group" not in findings_json["atomic_checks"]["functional_groups"]:
                        findings_json["atomic_checks"]["functional_groups"].append("Trifluoro group")
                    # Record the functional group introduction event
                    if "functional_group_introduction" not in findings_json["atomic_checks"]["named_reactions"]:
                        findings_json["atomic_checks"]["named_reactions"].append("functional_group_introduction")

            except Exception as e:
                # In a production setting, you might want to log this error
                pass

        for child in node.get("children", []):
            new_depth = depth
            if node["type"] != "reaction": # If current node is not a reaction (e.g., chemical), increase depth
                new_depth = depth + 1
            dfs_traverse(child, new_depth)

    dfs_traverse(route)

    # Check if CF3 is introduced early (in the first half of the synthesis)
    if cf3_introduction_depths:
        earliest_cf3_introduction = max(
            cf3_introduction_depths
        )  # Higher depth = earlier in synthesis

        # Early introduction means depth is greater than half of max_depth
        if earliest_cf3_introduction > max_depth / 2:
            result = True
            findings_json["structural_constraints"].append({
                "type": "positional",
                "details": {
                    "target_event": "functional_group_introduction",
                    "target_group": "Trifluoro group",
                    "position": "first_half"
                }
            })
        else:
            result = False
    elif target_has_cf3:
        # CF3 is in target but never introduced in any reaction - must be in starting materials
        result = True  # CF3 is present from the beginning, so it's introduced early
        findings_json["structural_constraints"].append({
            "type": "positional",
            "details": {
                "target_event": "functional_group_introduction",
                "target_group": "Trifluoro group",
                "position": "first_half"
            }
        })

    return result, findings_json
