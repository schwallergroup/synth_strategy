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
    Detects if an amine is introduced in the final step of the synthesis.
    """
    found_final_amine_introduction = False

    findings_template = {
        "atomic_checks": {
            "named_reactions": [],
            "ring_systems": [],
            "functional_groups": []
        },
        "structural_constraints": []
    }
    findings_json = copy.deepcopy(findings_template)

    print(f"Starting analysis of synthesis route")

    def dfs_traverse(node, depth=0):
        nonlocal found_final_amine_introduction, findings_json

        print(f"Traversing node at depth {depth}, type: {node['type']}")

        # Check if this is the final reaction (depth 1, child of the root product)
        if node["type"] == "reaction" and depth == 1:
            rsmi = node.get("metadata", {}).get("rsmi", "")
            if not rsmi:
                print("No reaction SMILES found")
                return

            # Extract reactants and product
            try:
                reactants_smiles = rsmi.split(">")[0].split(".")
                product_smiles = rsmi.split(">")[-1]

                print(f"Final reaction: {rsmi}")
                print(f"Reactants: {reactants_smiles}")
                print(f"Product: {product_smiles}")

                # Check if product contains amine groups
                has_primary_amine_in_product = checker.check_fg("Primary amine", product_smiles)
                if has_primary_amine_in_product:
                    findings_json["atomic_checks"]["functional_groups"].append("Primary amine")
                has_secondary_amine_in_product = checker.check_fg("Secondary amine", product_smiles)
                if has_secondary_amine_in_product:
                    findings_json["atomic_checks"]["functional_groups"].append("Secondary amine")
                has_tertiary_amine_in_product = checker.check_fg("Tertiary amine", product_smiles)
                if has_tertiary_amine_in_product:
                    findings_json["atomic_checks"]["functional_groups"].append("Tertiary amine")

                print(f"Primary amine in product: {has_primary_amine_in_product}")
                print(f"Secondary amine in product: {has_secondary_amine_in_product}")
                print(f"Tertiary amine in product: {has_tertiary_amine_in_product}")

                # Check if any reactant contains amine groups
                reactants_with_primary_amine = []
                for r in reactants_smiles:
                    if checker.check_fg("Primary amine", r):
                        reactants_with_primary_amine.append(True)
                        findings_json["atomic_checks"]["functional_groups"].append("Primary amine")
                    else:
                        reactants_with_primary_amine.append(False)

                reactants_with_secondary_amine = []
                for r in reactants_smiles:
                    if checker.check_fg("Secondary amine", r):
                        reactants_with_secondary_amine.append(True)
                        findings_json["atomic_checks"]["functional_groups"].append("Secondary amine")
                    else:
                        reactants_with_secondary_amine.append(False)

                reactants_with_tertiary_amine = []
                for r in reactants_smiles:
                    if checker.check_fg("Tertiary amine", r):
                        reactants_with_tertiary_amine.append(True)
                        findings_json["atomic_checks"]["functional_groups"].append("Tertiary amine")
                    else:
                        reactants_with_tertiary_amine.append(False)

                print(f"Reactants with primary amine: {reactants_with_primary_amine}")
                print(f"Reactants with secondary amine: {reactants_with_secondary_amine}")
                print(f"Reactants with tertiary amine: {reactants_with_tertiary_amine}")

                # Check for amine-introducing reactions
                is_reductive_amination = checker.check_reaction("{reductive amination}", rsmi)
                if is_reductive_amination:
                    findings_json["atomic_checks"]["named_reactions"].append("{reductive amination}")

                is_amide_reduction = (
                    checker.check_reaction("Reduction of primary amides to amines", rsmi)
                    or checker.check_reaction("Reduction of secondary amides to amines", rsmi)
                    or checker.check_reaction("Reduction of tertiary amides to amines", rsmi)
                )
                if checker.check_reaction("Reduction of primary amides to amines", rsmi):
                    findings_json["atomic_checks"]["named_reactions"].append("Reduction of primary amides to amines")
                if checker.check_reaction("Reduction of secondary amides to amines", rsmi):
                    findings_json["atomic_checks"]["named_reactions"].append("Reduction of secondary amides to amines")
                if checker.check_reaction("Reduction of tertiary amides to amines", rsmi):
                    findings_json["atomic_checks"]["named_reactions"].append("Reduction of tertiary amides to amines")

                is_nitrile_reduction = checker.check_reaction("Reduction of nitrile to amine", rsmi)
                if is_nitrile_reduction:
                    findings_json["atomic_checks"]["named_reactions"].append("Reduction of nitrile to amine")

                is_nitro_reduction = checker.check_reaction(
                    "Reduction of nitro groups to amines", rsmi
                )
                if is_nitro_reduction:
                    findings_json["atomic_checks"]["named_reactions"].append("Reduction of nitro groups to amines")

                is_azide_reduction = checker.check_reaction(
                    "Azide to amine reduction (Staudinger)", rsmi
                )
                if is_azide_reduction:
                    findings_json["atomic_checks"]["named_reactions"].append("Azide to amine reduction (Staudinger)")

                is_alkylation = checker.check_reaction("Alkylation of amines", rsmi)
                if is_alkylation:
                    findings_json["atomic_checks"]["named_reactions"].append("Alkylation of amines")

                is_buchwald_hartwig = checker.check_reaction("{Buchwald-Hartwig}", rsmi)
                if is_buchwald_hartwig:
                    findings_json["atomic_checks"]["named_reactions"].append("{Buchwald-Hartwig}")

                is_deprotection = (
                    checker.check_reaction("Boc amine deprotection", rsmi)
                    or checker.check_reaction("Tert-butyl deprotection of amine", rsmi)
                    or checker.check_reaction("Phthalimide deprotection", rsmi)
                    or checker.check_reaction("N-glutarimide deprotection", rsmi)
                )
                if checker.check_reaction("Boc amine deprotection", rsmi):
                    findings_json["atomic_checks"]["named_reactions"].append("Boc amine deprotection")
                if checker.check_reaction("Tert-butyl deprotection of amine", rsmi):
                    findings_json["atomic_checks"]["named_reactions"].append("Tert-butyl deprotection of amine")
                if checker.check_reaction("Phthalimide deprotection", rsmi):
                    findings_json["atomic_checks"]["named_reactions"].append("Phthalimide deprotection")
                if checker.check_reaction("N-glutarimide deprotection", rsmi):
                    findings_json["atomic_checks"]["named_reactions"].append("N-glutarimide deprotection")

                # Check if any amine-introducing reaction is detected
                amine_introducing_reaction = (
                    is_reductive_amination
                    or is_amide_reduction
                    or is_nitrile_reduction
                    or is_nitro_reduction
                    or is_azide_reduction
                    or is_alkylation
                    or is_buchwald_hartwig
                    or is_deprotection
                )

                print(f"Amine-introducing reaction detected: {amine_introducing_reaction}")
                print(
                    f"Reaction types: reductive amination: {is_reductive_amination}, amide reduction: {is_amide_reduction}, "
                    f"nitrile reduction: {is_nitrile_reduction}, nitro reduction: {is_nitro_reduction}, "
                    f"azide reduction: {is_azide_reduction}, alkylation: {is_alkylation}, "
                    f"Buchwald-Hartwig: {is_buchwald_hartwig}, deprotection: {is_deprotection}"
                )

                # Check if there are amines in product but not in reactants
                any_primary_amine_in_reactants = any(reactants_with_primary_amine)
                any_secondary_amine_in_reactants = any(reactants_with_secondary_amine)
                any_tertiary_amine_in_reactants = any(reactants_with_tertiary_amine)

                # Amine is introduced if:
                # 1. Product has an amine AND
                # 2. Either: a) Reactants don't have that type of amine OR b) It's an amine-introducing reaction
                primary_amine_introduced = has_primary_amine_in_product and (
                    not any_primary_amine_in_reactants or amine_introducing_reaction
                )
                secondary_amine_introduced = has_secondary_amine_in_product and (
                    not any_secondary_amine_in_reactants or amine_introducing_reaction
                )
                tertiary_amine_introduced = has_tertiary_amine_in_product and (
                    not any_tertiary_amine_in_reactants or amine_introducing_reaction
                )

                if (
                    primary_amine_introduced
                    or secondary_amine_introduced
                    or tertiary_amine_introduced
                ):
                    print(f"Amine introduction detected in final step")
                    print(
                        f"Primary amine in product: {has_primary_amine_in_product}, in reactants: {any_primary_amine_in_reactants}"
                    )
                    print(
                        f"Secondary amine in product: {has_secondary_amine_in_product}, in reactants: {any_secondary_amine_in_reactants}"
                    )
                    print(
                        f"Tertiary amine in product: {has_tertiary_amine_in_product}, in reactants: {any_tertiary_amine_in_reactants}"
                    )
                    print(f"Amine-introducing reaction detected: {amine_introducing_reaction}")
                    found_final_amine_introduction = True
                    # Add the structural constraint if the condition is met
                    findings_json["structural_constraints"].append({
                        "type": "positional",
                        "details": {
                            "target": "amine_introduction",
                            "position": "last_stage"
                        }
                    })

            except Exception as e:
                print(f"Error processing reaction: {e}")

        # Determine the new depth for the recursive call
        new_depth = depth
        if node['type'] != 'reaction': # Only increase depth when going from chemical to reaction
            new_depth = depth + 1

        # Process children with the determined new depth
        for child in node.get("children", []):
            dfs_traverse(child, new_depth)

    # Start traversal
    dfs_traverse(route)

    # Remove duplicates from atomic_checks lists
    for key in findings_json["atomic_checks"]:
        findings_json["atomic_checks"][key] = list(set(findings_json["atomic_checks"][key]))

    print(f"Final result: {found_final_amine_introduction}")
    return found_final_amine_introduction, findings_json
