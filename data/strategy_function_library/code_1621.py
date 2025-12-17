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


VALID_AMIDE_COUPLING_REACTIONS = [
    # General, robust checkers
    "Acylation of Nitrogen Nucleophiles by Carboxylic Acids",
    "Acylation of Nitrogen Nucleophiles by Acyl/Thioacyl/Carbamoyl Halides and Analogs_N",
    "Acylation of Nitrogen Nucleophiles by Acyl/Thioacyl/Carbamoyl Halides and Analogs_OS",
    "Carboxylic acid to amide conversion",
    # Specific, but correct and common named reactions
    "Schotten-Baumann_amide",
    # Specific checks for aminolysis of esters
    "Ester with primary amine to amide",
    "Ester with ammonia to amide",
    "Ester with secondary amine to amide"
]

def is_amide_formation_reaction(rsmi, reactants_smiles, product_smiles, findings_json: Dict) -> Tuple[bool, Dict]:
    """Helper function to check if a reaction is a valid amide formation.
    Also updates findings_json with detected atomic checks.
    """
    current_findings_json = copy.deepcopy(findings_json)
    
    has_amide_in_product = False
    if checker.check_fg("Primary amide", product_smiles):
        has_amide_in_product = True
        current_findings_json["atomic_checks"]["functional_groups"].append("Primary amide")
    if checker.check_fg("Secondary amide", product_smiles):
        has_amide_in_product = True
        current_findings_json["atomic_checks"]["functional_groups"].append("Secondary amide")
    if checker.check_fg("Tertiary amide", product_smiles):
        has_amide_in_product = True
        current_findings_json["atomic_checks"]["functional_groups"].append("Tertiary amide")

    if not has_amide_in_product:
        return False, current_findings_json

    # Check for amide in reactants to ensure it's newly formed
    has_amide_in_reactants = False
    for r in reactants_smiles:
        if checker.check_fg("Primary amide", r):
            has_amide_in_reactants = True
            current_findings_json["atomic_checks"]["functional_groups"].append("Primary amide")
        if checker.check_fg("Secondary amide", r):
            has_amide_in_reactants = True
            current_findings_json["atomic_checks"]["functional_groups"].append("Secondary amide")
        if checker.check_fg("Tertiary amide", r):
            has_amide_in_reactants = True
            current_findings_json["atomic_checks"]["functional_groups"].append("Tertiary amide")

    if has_amide_in_reactants:
        current_findings_json["structural_constraints"].append({"type": "negation", "details": {"target": "amide_functional_group", "scope": "reactants"}})
        return False, current_findings_json

    # Check if the reaction matches a known amide coupling type
    is_amide_coupling = False
    for reaction_name in VALID_AMIDE_COUPLING_REACTIONS:
        if checker.check_reaction(reaction_name, rsmi):
            is_amide_coupling = True
            current_findings_json["atomic_checks"]["named_reactions"].append(reaction_name)
            break

    if is_amide_coupling and has_amide_in_product:
        current_findings_json["structural_constraints"].append({"type": "co-occurrence", "details": {"targets": ["amide_coupling_reaction", "amide_in_product"]}})

    return is_amide_coupling, current_findings_json

def main(route) -> Tuple[bool, Dict]:
    """
    This function detects if the synthesis has a late-stage amide coupling
    (amide formation in the final step or penultimate step).
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

    has_late_stage_amide = False

    starts_with_molecule = route["type"] == "mol"

    def dfs_traverse(node, depth=0):
        nonlocal has_late_stage_amide, findings_json

        is_late_stage = (starts_with_molecule and depth <= 2) or (
            not starts_with_molecule and depth <= 1
        )

        if node.get("type") == "reaction" and is_late_stage:
            try:
                rsmi = node["metadata"]["mapped_reaction_smiles"]
                reactants_smiles = rsmi.split(">")[0].split(".")
                product_smiles = rsmi.split(">")[-1]

                amide_formation_result, updated_findings = is_amide_formation_reaction(rsmi, reactants_smiles, product_smiles, findings_json)
                findings_json = updated_findings # Update findings_json with results from helper

                if amide_formation_result:
                    has_late_stage_amide = True
                    # Add positional constraint if late stage and amide formation is detected
                    findings_json["structural_constraints"].append({"type": "positional", "details": {"target": "amide_formation", "position": "late_stage"}})
            except (KeyError, IndexError):
                pass

        for child in node.get("children", []):
            # New logic for depth calculation
            new_depth = depth + 1 if node.get("type") != "reaction" else depth
            dfs_traverse(child, new_depth)

    dfs_traverse(route)

    return has_late_stage_amide, findings_json
