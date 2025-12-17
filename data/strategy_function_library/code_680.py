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
    This function detects if a synthetic route employs a late-stage nitro reduction strategy.
    The strategy is flagged if a nitro reduction occurs in the final step (depth=1)
    and a nitro group was present in the product of at least one earlier step (depth>1).
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

    # Initialize tracking variables
    has_late_stage_nitro_reduction = False
    nitro_present_in_intermediates = False

    def dfs_traverse(node, depth=0):
        nonlocal has_late_stage_nitro_reduction, nitro_present_in_intermediates, findings_json

        if node["type"] == "reaction":
            # Extract reaction SMILES
            rsmi = node.get("metadata", {}).get("rsmi", "")
            if not rsmi:
                return

            product_smiles = rsmi.split(">")[-1]

            try:
                # Check if this is a nitro reduction reaction using the specific checker
                is_nitro_reduction = checker.check_reaction(
                    "Reduction of nitro groups to amines", rsmi
                )

                if is_nitro_reduction:
                    findings_json["atomic_checks"]["named_reactions"].append("Reduction of nitro groups to amines")
                    # Per the problem definition, depth=1 is the final step.
                    # We check for the reduction in the final step only.
                    if depth == 1:
                        has_late_stage_nitro_reduction = True
                        findings_json["structural_constraints"].append({
                            "type": "positional",
                            "details": {
                                "target": "Reduction of nitro groups to amines",
                                "position": "last_stage"
                            }
                        })
                else:
                    # If not a reduction, check for nitro groups in the product (intermediate compound)
                    # This confirms the nitro group was carried through the synthesis.
                    if depth > 1 and checker.check_fg("Nitro group", product_smiles):
                        nitro_present_in_intermediates = True
                        findings_json["atomic_checks"]["functional_groups"].append("Nitro group")

            except Exception as e:
                print(f"Error processing reaction: {e}")

        # Recursively process children
        for child in node.get("children", []):
            # New logic for depth calculation:
            # Depth increases only when traversing from a chemical node to a reaction node.
            # Depth remains the same when traversing from a reaction node to a chemical node.
            new_depth = depth if node["type"] == "reaction" else depth + 1
            dfs_traverse(child, new_depth)

    # Start traversal from the root
    dfs_traverse(route)

    # The strategy is present if we have a late-stage nitro reduction
    # and a nitro group was present in earlier intermediates.
    strategy_present = has_late_stage_nitro_reduction and nitro_present_in_intermediates

    if strategy_present:
        # Add the sequence constraint if both conditions are met
        findings_json["structural_constraints"].append({
            "type": "sequence",
            "details": {
                "first_event": "Nitro group",
                "second_event": "Reduction of nitro groups to amines"
            }
        })

    return strategy_present, findings_json
