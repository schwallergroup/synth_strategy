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


AMINE_FUNCTIONALIZATION_REACTIONS = [
    # Acylation reactions
    "Acylation of primary amines",
    "Acylation of secondary amines",
    "Acylation of Nitrogen Nucleophiles by Acyl/Thioacyl/Carbamoyl Halides and Analogs_N",
    "Acyl chloride with primary amine to amide (Schotten-Baumann)",
    "Acyl chloride with secondary amine to amide",
    "Carboxylic acid with primary amine to amide",
    "Ester with primary amine to amide",
    "Ester with secondary amine to amide",
    # Sulfonamide reactions
    "Sulfonamide synthesis (Schotten-Baumann) primary amine",
    "Sulfonamide synthesis (Schotten-Baumann) secondary amine",
    "sulfon_amide",
]

def main(route) -> Tuple[bool, Dict]:
    """
    Detects late-stage amine functionalization, defined as reactions occurring in the final two steps of a synthesis (depth <= 2). It identifies this strategy by first checking if the reaction matches a predefined list of amine functionalization reactions (see AMINE_FUNCTIONALIZATION_REACTIONS). If no named reaction is matched, it performs a fallback check for the presence of a primary amine, secondary amine, or aniline in the reactants and the formation of an amide, sulfonamide, urea, thiourea, or carbamate in the product.
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

    found_late_amine_functionalization = False

    def dfs_traverse(node, depth=0):
        nonlocal found_late_amine_functionalization, findings_json

        if node["type"] == "reaction" and depth <= 2:  # Check late-stage reactions (low depth)
            if "mapped_reaction_smiles" in node["metadata"]:
                rsmi = node["metadata"]["mapped_reaction_smiles"]
                reactants = rsmi.split(">")[0].split(".")
                product = rsmi.split(">")[-1]

                # Check if this is an amine functionalization reaction
                is_amine_functionalization = False
                for reaction_type in AMINE_FUNCTIONALIZATION_REACTIONS:
                    if checker.check_reaction(reaction_type, rsmi):
                        is_amine_functionalization = True
                        if reaction_type not in findings_json["atomic_checks"]["named_reactions"]:
                            findings_json["atomic_checks"]["named_reactions"].append(reaction_type)
                        break

                # If not a predefined reaction, check for functional group transformations
                if not is_amine_functionalization:
                    # Check for amine in reactants
                    has_amine_reactant = False
                    for reactant in reactants:
                        if checker.check_fg("Primary amine", reactant):
                            has_amine_reactant = True
                            if "Primary amine" not in findings_json["atomic_checks"]["functional_groups"]:
                                findings_json["atomic_checks"]["functional_groups"].append("Primary amine")
                        if checker.check_fg("Secondary amine", reactant):
                            has_amine_reactant = True
                            if "Secondary amine" not in findings_json["atomic_checks"]["functional_groups"]:
                                findings_json["atomic_checks"]["functional_groups"].append("Secondary amine")
                        if checker.check_fg("Aniline", reactant):
                            has_amine_reactant = True
                            if "Aniline" not in findings_json["atomic_checks"]["functional_groups"]:
                                findings_json["atomic_checks"]["functional_groups"].append("Aniline")
                        if has_amine_reactant: # Optimization: if any amine is found, no need to check others for this reactant
                            break

                    # Check for functionalized amine in product
                    has_functionalized_amine_product = False
                    if checker.check_fg("Primary amide", product):
                        has_functionalized_amine_product = True
                        if "Primary amide" not in findings_json["atomic_checks"]["functional_groups"]:
                            findings_json["atomic_checks"]["functional_groups"].append("Primary amide")
                    if checker.check_fg("Secondary amide", product):
                        has_functionalized_amine_product = True
                        if "Secondary amide" not in findings_json["atomic_checks"]["functional_groups"]:
                            findings_json["atomic_checks"]["functional_groups"].append("Secondary amide")
                    if checker.check_fg("Tertiary amide", product):
                        has_functionalized_amine_product = True
                        if "Tertiary amide" not in findings_json["atomic_checks"]["functional_groups"]:
                            findings_json["atomic_checks"]["functional_groups"].append("Tertiary amide")
                    if checker.check_fg("Sulfonamide", product):
                        has_functionalized_amine_product = True
                        if "Sulfonamide" not in findings_json["atomic_checks"]["functional_groups"]:
                            findings_json["atomic_checks"]["functional_groups"].append("Sulfonamide")
                    if checker.check_fg("Urea", product):
                        has_functionalized_amine_product = True
                        if "Urea" not in findings_json["atomic_checks"]["functional_groups"]:
                            findings_json["atomic_checks"]["functional_groups"].append("Urea")
                    if checker.check_fg("Thiourea", product):
                        has_functionalized_amine_product = True
                        if "Thiourea" not in findings_json["atomic_checks"]["functional_groups"]:
                            findings_json["atomic_checks"]["functional_groups"].append("Thiourea")
                    if checker.check_fg("Carbamate", product):
                        has_functionalized_amine_product = True
                        if "Carbamate" not in findings_json["atomic_checks"]["functional_groups"]:
                            findings_json["atomic_checks"]["functional_groups"].append("Carbamate")

                    # If we have both amine in reactants and functionalized amine in product
                    if has_amine_reactant and has_functionalized_amine_product:
                        is_amine_functionalization = True

                if is_amine_functionalization:
                    found_late_amine_functionalization = True
                    # Add structural constraint if late-stage amine functionalization is found
                    if {"type": "positional", "details": {"target": "amine_functionalization", "position": "late_stage (depth <= 2)"}} not in findings_json["structural_constraints"]:
                        findings_json["structural_constraints"].append({"type": "positional", "details": {"target": "amine_functionalization", "position": "late_stage (depth <= 2)"}})
                    return

        # Continue traversing
        for child in node.get("children", []):
            if not found_late_amine_functionalization:  # Stop traversal if already found
                # New depth calculation logic
                new_depth = depth
                if node["type"] != "reaction": # If current node is chemical, depth increases
                    new_depth = depth + 1
                # If current node is reaction, depth remains the same
                
                dfs_traverse(child, new_depth)

    # Start traversal from the root
    dfs_traverse(route)
    return found_late_amine_functionalization, findings_json
