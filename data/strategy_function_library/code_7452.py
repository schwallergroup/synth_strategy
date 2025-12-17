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


def main(route) -> Tuple[bool, Dict]:
    """
    This function detects if a trifluoromethyl group is preserved throughout the synthesis.
    It checks if the final product has a CF3 group and if this group is preserved
    in all steps after its introduction.
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

    final_product_has_cf3 = False

    # Track CF3 introduction and preservation along specific synthesis paths
    def dfs_traverse(node, depth=0, parent_has_cf3=False):
        nonlocal final_product_has_cf3, findings_json

        # Current node CF3 status
        current_has_cf3 = False

        if node["type"] == "mol":
            mol_smiles = node["smiles"]
            if checker.check_fg("Trifluoro group", mol_smiles):
                current_has_cf3 = True
                if "Trifluoro group" not in findings_json["atomic_checks"]["functional_groups"]:
                    findings_json["atomic_checks"]["functional_groups"].append("Trifluoro group")

            # Check if this is the final product (depth 0)
            if depth == 0:
                final_product_has_cf3 = current_has_cf3

            # Check for CF3 preservation only if parent had CF3
            elif parent_has_cf3 and not current_has_cf3:
                return False

        elif node["type"] == "reaction":
            # For reactions, we need to check if CF3 is preserved from product to reactants
            # (since we're traversing in retrosynthetic direction)
            if parent_has_cf3:
                try:
                    rsmi = node["metadata"]["mapped_reaction_smiles"]
                    reactants = rsmi.split(">")[0].split(".")
                    product = rsmi.split(">")[-1]

                    # In retrosynthesis, the product is the parent molecule (which has CF3)
                    # At least one reactant should have CF3 for preservation
                    reactant_has_cf3 = any(
                        checker.check_fg("Trifluoro group", r) for r in reactants
                    )
                    if reactant_has_cf3 and "Trifluoro group" not in findings_json["atomic_checks"]["functional_groups"]:
                        findings_json["atomic_checks"]["functional_groups"].append("Trifluoro group")

                    # Set current_has_cf3 to True if at least one reactant has CF3
                    current_has_cf3 = reactant_has_cf3

                except Exception as e:
                    return False

        # Determine the new depth for recursive calls
        new_depth = depth
        if node["type"] != "reaction":  # If current node is 'mol' (chemical), depth increases
            new_depth = depth + 1

        # Process children
        all_paths_preserve_cf3 = True

        for child in node.get("children", []):
            # Pass down the CF3 status to children
            path_preserves_cf3 = dfs_traverse(child, new_depth, current_has_cf3)

            # Only evaluate preservation for paths where parent had CF3
            if current_has_cf3 and not path_preserves_cf3:
                all_paths_preserve_cf3 = False

        # Only check preservation for paths that actually contain CF3
        if parent_has_cf3:
            return all_paths_preserve_cf3

        # If no CF3 in parent path, return True (no preservation needed)
        return True

    # Start traversal from the root
    preservation_result = dfs_traverse(route)

    # Final result: product has CF3 and it was preserved throughout synthesis
    result = final_product_has_cf3 and preservation_result

    if final_product_has_cf3:
        findings_json["structural_constraints"].append({
            "type": "positional",
            "details": {
                "target": "Trifluoro group",
                "position": "last_stage"
            }
        })

    if preservation_result:
        # This implies no loss of Trifluoro group
        findings_json["structural_constraints"].append({
            "type": "negation",
            "details": {
                "target": "loss_of_Trifluoro group"
            }
        })

    return result, findings_json