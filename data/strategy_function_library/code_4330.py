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


SULFONATE_LEAVING_GROUPS = ["Mesylate", "Triflate", "Tosylate"]
ALCOHOL_TYPES = ["Primary alcohol", "Secondary alcohol", "Tertiary alcohol", "Aromatic alcohol"]

def main(route) -> Tuple[bool, Dict]:
    """
    Detects a three-step sequence for introducing an azide: 1) formation of an alcohol, 2) activation of the alcohol into a sulfonate ester, and 3) substitution with an azide nucleophile. The specific sulfonate esters checked are defined in SULFONATE_LEAVING_GROUPS, and the alcohol types are defined in ALCOHOL_TYPES.
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

    # Track if we found each step in the sequence
    found_azide_introduction = False
    found_alcohol_activation = False
    found_alcohol_formation = False

    # Track reaction sequence for order verification
    reaction_sequence = []

    def dfs_traverse(node, depth=0):
        nonlocal found_azide_introduction, found_alcohol_activation, found_alcohol_formation, findings_json

        if node["type"] == "reaction":
            # Extract reactants and product
            rsmi = node["metadata"]["mapped_reaction_smiles"]
            reactants_smiles = rsmi.split(">")[0].split(".")
            product_smiles = rsmi.split(">")[-1]

            # Check for azide introduction
            product_has_azide = checker.check_fg("Azide", product_smiles)
            if product_has_azide:
                if "Azide" not in findings_json["atomic_checks"]["functional_groups"]:
                    findings_json["atomic_checks"]["functional_groups"].append("Azide")

            reactants_have_azide = any(checker.check_fg("Azide", r) for r in reactants_smiles if r)
            reactants_have_leaving_group = False
            for lg in SULFONATE_LEAVING_GROUPS:
                if any(checker.check_fg(lg, r) for r in reactants_smiles if r):
                    reactants_have_leaving_group = True
                    if lg not in findings_json["atomic_checks"]["functional_groups"]:
                        findings_json["atomic_checks"]["functional_groups"].append(lg)

            if product_has_azide and not reactants_have_azide and reactants_have_leaving_group:
                found_azide_introduction = True
                reaction_sequence.append(("azide_introduction", depth))
                if "azide_introduction" not in findings_json["atomic_checks"]["named_reactions"]:
                    findings_json["atomic_checks"]["named_reactions"].append("azide_introduction")

            # Check for sulfonate activation
            product_has_leaving_group = False
            for lg in SULFONATE_LEAVING_GROUPS:
                if checker.check_fg(lg, product_smiles):
                    product_has_leaving_group = True
                    if lg not in findings_json["atomic_checks"]["functional_groups"]:
                        findings_json["atomic_checks"]["functional_groups"].append(lg)

            reactants_have_leaving_group_for_activation = False
            for lg in SULFONATE_LEAVING_GROUPS:
                if any(checker.check_fg(lg, r) for r in reactants_smiles if r):
                    reactants_have_leaving_group_for_activation = True
                    if lg not in findings_json["atomic_checks"]["functional_groups"]:
                        findings_json["atomic_checks"]["functional_groups"].append(lg)

            reactants_have_alcohol = False
            for alc in ALCOHOL_TYPES:
                if any(checker.check_fg(alc, r) for r in reactants_smiles if r):
                    reactants_have_alcohol = True
                    if alc not in findings_json["atomic_checks"]["functional_groups"]:
                        findings_json["atomic_checks"]["functional_groups"].append(alc)

            if (
                product_has_leaving_group
                and not reactants_have_leaving_group_for_activation
                and reactants_have_alcohol
            ):
                found_alcohol_activation = True
                reaction_sequence.append(("alcohol_activation", depth))
                if "alcohol_activation" not in findings_json["atomic_checks"]["named_reactions"]:
                    findings_json["atomic_checks"]["named_reactions"].append("alcohol_activation")

            # Check for alcohol formation
            product_has_alcohol = False
            for alc in ALCOHOL_TYPES:
                if checker.check_fg(alc, product_smiles):
                    product_has_alcohol = True
                    if alc not in findings_json["atomic_checks"]["functional_groups"]:
                        findings_json["atomic_checks"]["functional_groups"].append(alc)

            reactants_have_alcohol_for_formation = False
            for alc in ALCOHOL_TYPES:
                if any(checker.check_fg(alc, r) for r in reactants_smiles if r):
                    reactants_have_alcohol_for_formation = True
                    if alc not in findings_json["atomic_checks"]["functional_groups"]:
                        findings_json["atomic_checks"]["functional_groups"].append(alc)

            if product_has_alcohol and not reactants_have_alcohol_for_formation:
                found_alcohol_formation = True
                reaction_sequence.append(("alcohol_formation", depth))
                if "alcohol_formation" not in findings_json["atomic_checks"]["named_reactions"]:
                    findings_json["atomic_checks"]["named_reactions"].append("alcohol_formation")

        # Continue traversing
        for child in node.get("children", []):
            # New logic for depth calculation
            new_depth = depth if node["type"] == "reaction" else depth + 1
            dfs_traverse(child, new_depth)

    # Start traversal
    dfs_traverse(route)

    # Check if all steps in the strategy are found
    all_steps_found = (
        found_azide_introduction and found_alcohol_activation and found_alcohol_formation
    )

    # Check sequence order (in retrosynthetic direction: azide → leaving group → alcohol)
    correct_sequence = False
    if reaction_sequence:
        # Sort by depth (higher depth = earlier in synthesis)
        reaction_sequence.sort(key=lambda x: x[1], reverse=True)
        reaction_types = [r[0] for r in reaction_sequence]

        # Check if the sequence contains the expected pattern
        for i in range(len(reaction_types) - 2):
            if (
                reaction_types[i] == "alcohol_formation"
                and reaction_types[i + 1] == "alcohol_activation"
                and reaction_types[i + 2] == "azide_introduction"
            ):
                correct_sequence = True
                break

    result = all_steps_found and correct_sequence

    if all_steps_found:
        # Add co-occurrence constraint if all steps are found
        findings_json["structural_constraints"].append({
            "type": "co-occurrence",
            "details": {
                "targets": [
                    "azide_introduction",
                    "alcohol_activation",
                    "alcohol_formation"
                ]
            }
        })

    if correct_sequence:
        # Add sequence constraint if the correct order is found
        findings_json["structural_constraints"].append({
            "type": "sequence",
            "details": {
                "ordered_events": [
                    "alcohol_formation",
                    "alcohol_activation",
                    "azide_introduction"
                ],
                "direction": "forward_synthesis"
            }
        })

    return result, findings_json
