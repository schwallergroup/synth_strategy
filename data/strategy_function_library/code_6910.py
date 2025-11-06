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
    Detects a potential functional group relay sequence (Alcohol → Mesylate → Azide → Amine) within a synthesis route.
    It verifies that each transformation type is present and that they occur in the correct synthetic order based on their depth in the route.
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

    # Track the transformations and their sequence
    transformations = []

    def dfs_traverse(node, depth=0):
        nonlocal transformations, findings_json
        if node["type"] == "reaction":
            # Extract reactants and product
            rsmi = node["metadata"]["mapped_reaction_smiles"]
            reactants_smiles = rsmi.split(">")[0].split(".")
            product_smiles = rsmi.split(">")[-1]

            # Check for alcohol to mesylate transformation
            alcohol_found_in_reactants = False
            for r in reactants_smiles:
                if checker.check_fg("Primary alcohol", r):
                    findings_json["atomic_checks"]["functional_groups"].append("Primary alcohol")
                    alcohol_found_in_reactants = True
                if checker.check_fg("Secondary alcohol", r):
                    findings_json["atomic_checks"]["functional_groups"].append("Secondary alcohol")
                    alcohol_found_in_reactants = True
                if checker.check_fg("Tertiary alcohol", r):
                    findings_json["atomic_checks"]["functional_groups"].append("Tertiary alcohol")
                    alcohol_found_in_reactants = True
            mesylate_found_in_product = checker.check_fg("Mesylate", product_smiles)
            if mesylate_found_in_product:
                findings_json["atomic_checks"]["functional_groups"].append("Mesylate")

            if alcohol_found_in_reactants and mesylate_found_in_product:
                transformations.append(("alcohol_to_mesylate", depth))
                findings_json["atomic_checks"]["named_reactions"].append("alcohol_to_mesylate")

            # Check for mesylate to azide transformation
            mesylate_found_in_reactants = False
            for r in reactants_smiles:
                if checker.check_fg("Mesylate", r):
                    findings_json["atomic_checks"]["functional_groups"].append("Mesylate")
                    mesylate_found_in_reactants = True
            azide_found_in_product = checker.check_fg("Azide", product_smiles)
            if azide_found_in_product:
                findings_json["atomic_checks"]["functional_groups"].append("Azide")

            if mesylate_found_in_reactants and azide_found_in_product:
                transformations.append(("mesylate_to_azide", depth))
                findings_json["atomic_checks"]["named_reactions"].append("mesylate_to_azide")

            # Check for azide to amine transformation
            azide_found_in_reactants = False
            for r in reactants_smiles:
                if checker.check_fg("Azide", r):
                    findings_json["atomic_checks"]["functional_groups"].append("Azide")
                    azide_found_in_reactants = True
            amine_found_in_product = False
            if checker.check_fg("Primary amine", product_smiles):
                findings_json["atomic_checks"]["functional_groups"].append("Primary amine")
                amine_found_in_product = True
            if checker.check_fg("Secondary amine", product_smiles):
                findings_json["atomic_checks"]["functional_groups"].append("Secondary amine")
                amine_found_in_product = True
            if checker.check_fg("Tertiary amine", product_smiles):
                findings_json["atomic_checks"]["functional_groups"].append("Tertiary amine")
                amine_found_in_product = True

            if azide_found_in_reactants and amine_found_in_product:
                transformations.append(("azide_to_amine", depth))
                findings_json["atomic_checks"]["named_reactions"].append("azide_to_amine")

        # Traverse children
        for child in node.get("children", []):
            # New logic for depth calculation
            new_depth = depth
            if node["type"] != "reaction":  # If current node is chemical, depth increases
                new_depth = depth + 1
            dfs_traverse(child, new_depth)

    # Start traversal
    dfs_traverse(route)

    # Check if all transformations are present and in the correct sequence
    has_alcohol_to_mesylate = any(t[0] == "alcohol_to_mesylate" for t in transformations)
    has_mesylate_to_azide = any(t[0] == "mesylate_to_azide" for t in transformations)
    has_azide_to_amine = any(t[0] == "azide_to_amine" for t in transformations)

    # Check sequence by comparing depths
    in_correct_sequence = False
    if has_alcohol_to_mesylate and has_mesylate_to_azide and has_azide_to_amine:
        alcohol_to_mesylate_depth = next(
            t[1] for t in transformations if t[0] == "alcohol_to_mesylate"
        )
        mesylate_to_azide_depth = next(t[1] for t in transformations if t[0] == "mesylate_to_azide")
        azide_to_amine_depth = next(t[1] for t in transformations if t[0] == "azide_to_amine")

        # In retrosynthetic traversal, later steps have lower depth
        in_correct_sequence = (
            azide_to_amine_depth < mesylate_to_azide_depth < alcohol_to_mesylate_depth
        )

    relay_strategy_present = (
        has_alcohol_to_mesylate
        and has_mesylate_to_azide
        and has_azide_to_amine
        and in_correct_sequence
    )

    if relay_strategy_present:
        # Add co-occurrence constraint
        findings_json["structural_constraints"].append({
            "type": "co-occurrence",
            "details": {
                "targets": [
                    "alcohol_to_mesylate",
                    "mesylate_to_azide",
                    "azide_to_amine"
                ]
            }
        })
        # Add sequence constraint
        findings_json["structural_constraints"].append({
            "type": "sequence",
            "details": {
                "ordered_events": [
                    "alcohol_to_mesylate",
                    "mesylate_to_azide",
                    "azide_to_amine"
                ],
                "description": "The sequence of transformations must follow Alcohol -> Mesylate, then Mesylate -> Azide, then Azide -> Amine in the forward synthesis direction."
            }
        })

    return relay_strategy_present, findings_json