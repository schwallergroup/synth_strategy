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


AROMATIC_NITRATION_REACTIONS = [
    "Aromatic nitration with HNO3",
    "Aromatic nitration with NO3 salt",
    "Aromatic nitration with NO2 salt",
    "Aromatic nitration with alkyl NO2",
]

UREA_FORMATION_REACTIONS = [
    "Urea synthesis via isocyanate and primary amine",
    "Urea synthesis via isocyanate and secondary amine",
    "Urea synthesis via isocyanate and diazo",
    "Urea synthesis via isocyanate and sulfonamide",
    "{urea}",
]

def main(route) -> Tuple[bool, Dict]:
    """
    Detects a multi-step synthetic strategy involving the sequence: 1) Aromatic Nitration, 2) Nitro Group Reduction, and 3) Urea Formation. This function verifies that these three reaction types occur in the correct synthetic order (nitration must precede reduction, which must precede urea formation). It identifies specific, named reaction types for nitration and urea formation by checking against the `AROMATIC_NITRATION_REACTIONS` and `UREA_FORMATION_REACTIONS` lists, respectively.
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
    has_nitration = False
    has_nitro_reduction = False
    has_urea_formation = False
    nitration_depth = -1
    nitro_reduction_depth = -1
    urea_formation_depth = -1

    def dfs_traverse(node, depth=0):
        nonlocal has_nitration, has_nitro_reduction, has_urea_formation
        nonlocal nitration_depth, nitro_reduction_depth, urea_formation_depth
        nonlocal findings_json

        if node["type"] == "reaction":
            # Extract reaction SMILES
            rsmi = node["metadata"].get("rsmi", "")
            if not rsmi:
                return

            # Check for nitration reaction
            for r_name in AROMATIC_NITRATION_REACTIONS:
                if checker.check_reaction(r_name, rsmi):
                    has_nitration = True
                    nitration_depth = depth
                    if r_name not in findings_json["atomic_checks"]["named_reactions"]:
                        findings_json["atomic_checks"]["named_reactions"].append(r_name)
                    break

            # Check for nitro reduction
            if checker.check_reaction("Reduction of nitro groups to amines", rsmi):
                has_nitro_reduction = True
                nitro_reduction_depth = depth
                if "Reduction of nitro groups to amines" not in findings_json["atomic_checks"]["named_reactions"]:
                    findings_json["atomic_checks"]["named_reactions"].append("Reduction of nitro groups to amines")

            # Check for urea formation
            for r_name in UREA_FORMATION_REACTIONS:
                if checker.check_reaction(r_name, rsmi):
                    has_urea_formation = True
                    urea_formation_depth = depth
                    if r_name not in findings_json["atomic_checks"]["named_reactions"]:
                        findings_json["atomic_checks"]["named_reactions"].append(r_name)
                    break

        # Traverse children
        for child in node.get("children", []):
            # New logic: depth increases only when traversing from chemical to reaction node
            # Depth remains the same when traversing from reaction to chemical node
            if node["type"] == "reaction":
                dfs_traverse(child, depth)
            else:
                dfs_traverse(child, depth + 1)

    # Start traversal
    dfs_traverse(route)

    # Check if the strategy is present
    # In retrosynthetic traversal, higher depth means earlier in the synthetic route
    # So nitro reduction should have a higher depth than urea formation
    strategy_present = False

    if has_nitro_reduction and has_urea_formation and nitro_reduction_depth > urea_formation_depth:
        strategy_present = True
        # Add structural constraint for co-occurrence and sequence
        findings_json["structural_constraints"].append({
            "type": "co-occurrence",
            "details": {
                "targets": [
                    "Reduction of nitro groups to amines",
                    "Urea formation"
                ]
            }
        })
        findings_json["structural_constraints"].append({
            "type": "sequence",
            "details": {
                "before": "Reduction of nitro groups to amines",
                "after": "Urea formation"
            }
        })

    # Ideally, we should also have nitration before nitro reduction
    if has_nitration and strategy_present: # Only check if the main strategy is already present
        if nitration_depth > nitro_reduction_depth:
            strategy_present = True # This condition is met, so strategy remains true
            findings_json["structural_constraints"].append({
                "type": "sequence",
                "details": {
                    "before": "Aromatic nitration",
                    "after": "Reduction of nitro groups to amines"
                }
            })
        else:
            strategy_present = False # If nitration is present but order is wrong, strategy is false

    return strategy_present, findings_json
