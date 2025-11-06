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


ESTERIFICATION_REACTIONS = [
    "Esterification of Carboxylic Acids",
    "Schotten-Baumann to ester",
    "O-alkylation of carboxylic acids with diazo compounds",
    "Transesterification",
    "Oxidative esterification of primary alcohols",
]

def main(route) -> Tuple[bool, Dict]:
    """
    This function detects if the synthesis involves an esterification step.
    It identifies the reaction by checking for the formation of an ester from a carboxylic acid derivative and an alcohol,
    or by identifying a specific named reaction. The named reactions checked are:
    Esterification of Carboxylic Acids, Schotten-Baumann to ester, O-alkylation of carboxylic acids with diazo compounds, Transesterification, and Oxidative esterification of primary alcohols.
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

    has_esterification = False

    def dfs_traverse(node, current_depth=0):
        nonlocal has_esterification, findings_json

        if node["type"] == "reaction":
            rsmi = node.get("metadata", {}).get("rsmi", "")
            if rsmi:
                reactants = rsmi.split(">")[0].split(".")
                product = rsmi.split(">")[-1]

                # Check for carboxylic acid or activated derivatives in reactants
                has_acid = False
                for reactant in reactants:
                    if checker.check_fg("Carboxylic acid", reactant):
                        findings_json["atomic_checks"]["functional_groups"].append("Carboxylic acid")
                        has_acid = True
                    if checker.check_fg("Acyl halide", reactant):
                        findings_json["atomic_checks"]["functional_groups"].append("Acyl halide")
                        has_acid = True
                    if checker.check_fg("Anhydride", reactant):
                        findings_json["atomic_checks"]["functional_groups"].append("Anhydride")
                        has_acid = True

                # Check for alcohol in reactants
                has_alcohol = False
                for reactant in reactants:
                    if checker.check_fg("Primary alcohol", reactant):
                        findings_json["atomic_checks"]["functional_groups"].append("Primary alcohol")
                        has_alcohol = True
                    if checker.check_fg("Secondary alcohol", reactant):
                        findings_json["atomic_checks"]["functional_groups"].append("Secondary alcohol")
                        has_alcohol = True
                    if checker.check_fg("Tertiary alcohol", reactant):
                        findings_json["atomic_checks"]["functional_groups"].append("Tertiary alcohol")
                        has_alcohol = True
                    if checker.check_fg("Aromatic alcohol", reactant):
                        findings_json["atomic_checks"]["functional_groups"].append("Aromatic alcohol")
                        has_alcohol = True
                    if checker.check_fg("Phenol", reactant):
                        findings_json["atomic_checks"]["functional_groups"].append("Phenol")
                        has_alcohol = True

                # Check for ester in product
                has_ester = False
                if checker.check_fg("Ester", product):
                    findings_json["atomic_checks"]["functional_groups"].append("Ester")
                    has_ester = True

                # Verify this is an esterification reaction (multiple pathways)
                is_esterification = False
                for rxn_name in ESTERIFICATION_REACTIONS:
                    if checker.check_reaction(rxn_name, rsmi):
                        findings_json["atomic_checks"]["named_reactions"].append(rxn_name)
                        is_esterification = True
                        break

                # Direct pattern check: if we have carboxylic acid/derivative + alcohol → ester
                if has_acid and has_alcohol and has_ester:
                    # Check if the ester in the product is new (not present in reactants)
                    ester_in_reactants = any(
                        checker.check_fg("Ester", reactant) for reactant in reactants
                    )
                    if not ester_in_reactants or is_esterification:
                        has_esterification = True
                        # This implies a structural constraint: formation of ester from acid/alcohol
                        # The original JSON has no specific structural constraint for this, so we add a generic one if needed.
                        # For this problem, we assume the overall 'has_esterification' flag is the primary structural constraint.
                        if {"type": "esterification_event", "details": {"description": "Ester formed from acid/derivative and alcohol"}} not in findings_json["structural_constraints"]:
                            findings_json["structural_constraints"].append({"type": "esterification_event", "details": {"description": "Ester formed from acid/derivative and alcohol"}})

                if is_esterification and not has_esterification: # If it's a named esterification but the acid/alcohol/ester pattern wasn't met (e.g., transesterification)
                    has_esterification = True
                    if {"type": "esterification_event", "details": {"description": "Named esterification reaction detected"}} not in findings_json["structural_constraints"]:
                        findings_json["structural_constraints"].append({"type": "esterification_event", "details": {"description": "Named esterification reaction detected"}})

        # Determine the next depth based on the current node's type
        next_depth = current_depth
        if node["type"] != "reaction":  # If current node is 'chemical' or other non-reaction type
            next_depth = current_depth + 1

        # Traverse children with the determined depth
        for child in node.get("children", []):
            dfs_traverse(child, next_depth)

    # Start traversal from the root
    dfs_traverse(route)

    # Remove duplicates from atomic_checks lists
    for key in findings_json["atomic_checks"]:
        findings_json["atomic_checks"][key] = list(set(findings_json["atomic_checks"][key]))

    return has_esterification, findings_json