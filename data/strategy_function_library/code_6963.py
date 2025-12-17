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


AMIDE_COUPLING_REACTIONS_FROM_DERIVATIVES = [
    "Acyl chloride with primary amine to amide (Schotten-Baumann)",
    "Ester with primary amine to amide",
    "Ester with secondary amine to amide",
    "Acyl chloride with secondary amine to amide",
]

def main(route) -> Tuple[bool, Dict]:
    """
    This function detects a synthetic strategy involving late-stage amide coupling
    (final reaction forms an amide bond between carboxylic acid and amine).
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

    amide_coupling_at_depth_0 = False

    def dfs_traverse(node, depth=0):
        nonlocal amide_coupling_at_depth_0, findings_json

        if node["type"] == "reaction" and "metadata" in node and "mapped_reaction_smiles" in node["metadata"]:
            rsmi = node["metadata"]["mapped_reaction_smiles"]

            # Extract depth information
            # The depth is now passed as an argument to the function

            print(f"Examining reaction at depth {depth}: {rsmi}")

            if depth == 1:  # Final reaction
                reactants = rsmi.split(">")[0].split(".")
                product = rsmi.split(">")[-1]

                print(f"Final reaction - Reactants: {reactants}, Product: {product}")

                # Check for carboxylic acid in reactants
                acid_present = False
                for r in reactants:
                    if checker.check_fg("Carboxylic acid", r):
                        acid_present = True
                        if "Carboxylic acid" not in findings_json["atomic_checks"]["functional_groups"]:
                            findings_json["atomic_checks"]["functional_groups"].append("Carboxylic acid")

                # Check for amine in reactants (primary or secondary)
                amine_present = False
                for r in reactants:
                    if checker.check_fg("Primary amine", r):
                        amine_present = True
                        if "Primary amine" not in findings_json["atomic_checks"]["functional_groups"]:
                            findings_json["atomic_checks"]["functional_groups"].append("Primary amine")
                    if checker.check_fg("Secondary amine", r):
                        amine_present = True
                        if "Secondary amine" not in findings_json["atomic_checks"]["functional_groups"]:
                            findings_json["atomic_checks"]["functional_groups"].append("Secondary amine")

                # Check for amide in product
                amide_in_product = False
                if checker.check_fg("Primary amide", product):
                    amide_in_product = True
                    if "Primary amide" not in findings_json["atomic_checks"]["functional_groups"]:
                        findings_json["atomic_checks"]["functional_groups"].append("Primary amide")
                if checker.check_fg("Secondary amide", product):
                    amide_in_product = True
                    if "Secondary amide" not in findings_json["atomic_checks"]["functional_groups"]:
                        findings_json["atomic_checks"]["functional_groups"].append("Secondary amide")
                if checker.check_fg("Tertiary amide", product):
                    amide_in_product = True
                    if "Tertiary amide" not in findings_json["atomic_checks"]["functional_groups"]:
                        findings_json["atomic_checks"]["functional_groups"].append("Tertiary amide")

                print(
                    f"Functional group checks - Acid: {acid_present}, Amine: {amine_present}, Amide in product: {amide_in_product}"
                )

                # Check for known amide coupling reactions
                is_amide_coupling = False
                for r_name in AMIDE_COUPLING_REACTIONS_FROM_DERIVATIVES:
                    if checker.check_reaction(r_name, rsmi):
                        is_amide_coupling = True
                        if r_name not in findings_json["atomic_checks"]["named_reactions"]:
                            findings_json["atomic_checks"]["named_reactions"].append(r_name)

                print(f"Is known amide coupling reaction: {is_amide_coupling}")

                if (acid_present and amine_present and amide_in_product):
                    amide_coupling_at_depth_0 = True
                    print(f"Detected generic amide coupling at final step (depth 1): {rsmi}")
                    # Add the generic amide coupling structural constraint
                    if {"type": "co-occurrence", "details": {"scope": "single_reaction_step", "position": "last_stage", "targets": ["reactant_has_fg:Carboxylic acid", "reactant_has_fg:Amine", "product_has_fg:Amide"], "description": "A generic amide coupling (acid + amine -> amide) occurs in the final reaction step."}} not in findings_json["structural_constraints"]:
                        findings_json["structural_constraints"].append({"type": "co-occurrence", "details": {"scope": "single_reaction_step", "position": "last_stage", "targets": ["reactant_has_fg:Carboxylic acid", "reactant_has_fg:Amine", "product_has_fg:Amide"], "description": "A generic amide coupling (acid + amine -> amide) occurs in the final reaction step."}})

                if is_amide_coupling:
                    amide_coupling_at_depth_0 = True
                    print(f"Detected specific amide coupling reaction at final step (depth 1): {rsmi}")
                    # Add the specific amide coupling structural constraint
                    if {"type": "positional", "details": {"position": "last_stage", "target_options": ["Acyl chloride with primary amine to amide (Schotten-Baumann)", "Ester with primary amine to amide", "Ester with secondary amine to amide", "Acyl chloride with secondary amine to amide"], "description": "The final reaction step is one of several known amide coupling reactions."}} not in findings_json["structural_constraints"]:
                        findings_json["structural_constraints"].append({"type": "positional", "details": {"position": "last_stage", "target_options": ["Acyl chloride with primary amine to amide (Schotten-Baumann)", "Ester with primary amine to amide", "Ester with secondary amine to amide", "Acyl chloride with secondary amine to amide"], "description": "The final reaction step is one of several known amide coupling reactions."}})

        for child in node.get("children", []):
            # New depth calculation logic
            new_depth = depth
            if node["type"] != "reaction": # If current node is chemical, depth increases
                new_depth = depth + 1
            # If current node is reaction, depth remains the same
            dfs_traverse(child, new_depth)

    dfs_traverse(route)
    print(f"Final result: {amide_coupling_at_depth_0}")
    return amide_coupling_at_depth_0, findings_json