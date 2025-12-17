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
    Detects if a synthetic route contains both an alcohol oxidation step (to a carbonyl compound like a ketone, aldehyde, or carboxylic acid) and a carbonyl reduction step (from a ketone or aldehyde to an alcohol). This identifies routes with bidirectional redox transformations of the alcohol/carbonyl functional group.
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

    alcohol_to_ketone = False
    ketone_to_alcohol = False

    print("Starting bidirectional alcohol-ketone interconversion check")

    def dfs_traverse(node, depth=0):
        nonlocal alcohol_to_ketone, ketone_to_alcohol, findings_json

        if node["type"] == "reaction":
            if "mapped_reaction_smiles" in node.get("metadata", {}):
                rsmi = node["metadata"]["mapped_reaction_smiles"]
                reactants_str = rsmi.split(">")[0]
                product_str = rsmi.split(">")[-1]

                print(f"Checking reaction: {rsmi}")

                try:
                    # Split reactants for individual checking
                    reactants = reactants_str.split(".")

                    # Check for alcohol to ketone/aldehyde/carboxylic acid/ester (oxidation)
                    oxidation_reactions = [
                        "Oxidation or Dehydrogenation of Alcohols to Aldehydes and Ketones",
                        "Oxidation of alcohol to carboxylic acid",
                        "Oxidation of alcohol and aldehyde to ester",
                        "Oxidative esterification of primary alcohols"
                    ]
                    
                    current_oxidation_reaction = None
                    for r_name in oxidation_reactions:
                        if checker.check_reaction(r_name, rsmi):
                            current_oxidation_reaction = r_name
                            break

                    if current_oxidation_reaction:
                        # Verify that an alcohol is present in reactants
                        alcohol_fgs = [
                            "Primary alcohol", "Secondary alcohol", "Tertiary alcohol",
                            "Aromatic alcohol", "Enol"
                        ]
                        has_alcohol = False
                        for r in reactants:
                            for fg_name in alcohol_fgs:
                                if checker.check_fg(fg_name, r):
                                    if fg_name not in findings_json["atomic_checks"]["functional_groups"]:
                                        findings_json["atomic_checks"]["functional_groups"].append(fg_name)
                                    has_alcohol = True

                        # Verify that a carbonyl compound is present in product
                        carbonyl_fgs_product = [
                            "Ketone", "Aldehyde", "Carboxylic acid", "Ester"
                        ]
                        has_carbonyl = False
                        for fg_name in carbonyl_fgs_product:
                            if checker.check_fg(fg_name, product_str):
                                if fg_name not in findings_json["atomic_checks"]["functional_groups"]:
                                    findings_json["atomic_checks"]["functional_groups"].append(fg_name)
                                has_carbonyl = True

                        if has_alcohol and has_carbonyl:
                            print("Alcohol to carbonyl transformation detected")
                            alcohol_to_ketone = True
                            if current_oxidation_reaction not in findings_json["atomic_checks"]["named_reactions"]:
                                findings_json["atomic_checks"]["named_reactions"].append(current_oxidation_reaction)

                    # Check for ketone/aldehyde/ester to alcohol (reduction)
                    reduction_reactions = [
                        "Reduction of aldehydes and ketones to alcohols",
                        "Reduction of ketone to secondary alcohol",
                        "Grignard from aldehyde to alcohol",
                        "Grignard from ketone to alcohol"
                    ]

                    current_reduction_reaction = None
                    for r_name in reduction_reactions:
                        if checker.check_reaction(r_name, rsmi):
                            current_reduction_reaction = r_name
                            break

                    if current_reduction_reaction:
                        # Verify that a carbonyl compound is present in reactants
                        carbonyl_fgs_reactants = [
                            "Ketone", "Aldehyde", "Formaldehyde", "Ester", "Carboxylic acid"
                        ]
                        has_carbonyl = False
                        for r in reactants:
                            for fg_name in carbonyl_fgs_reactants:
                                if checker.check_fg(fg_name, r):
                                    if fg_name not in findings_json["atomic_checks"]["functional_groups"]:
                                        findings_json["atomic_checks"]["functional_groups"].append(fg_name)
                                    has_carbonyl = True

                        # Verify that an alcohol is present in product
                        alcohol_fgs_product = [
                            "Primary alcohol", "Secondary alcohol", "Tertiary alcohol",
                            "Aromatic alcohol"
                        ]
                        has_alcohol = False
                        for fg_name in alcohol_fgs_product:
                            if checker.check_fg(fg_name, product_str):
                                if fg_name not in findings_json["atomic_checks"]["functional_groups"]:
                                    findings_json["atomic_checks"]["functional_groups"].append(fg_name)
                                has_alcohol = True

                        if has_carbonyl and has_alcohol:
                            print("Carbonyl to alcohol transformation detected")
                            ketone_to_alcohol = True
                            if current_reduction_reaction not in findings_json["atomic_checks"]["named_reactions"]:
                                findings_json["atomic_checks"]["named_reactions"].append(current_reduction_reaction)

                except Exception as e:
                    print(f"Error processing reaction: {e}")

                print(
                    f"Current status - alcohol_to_ketone: {alcohol_to_ketone}, ketone_to_alcohol: {ketone_to_alcohol}"
                )

        # Continue traversal
        for child in node.get("children", []):
            new_depth = depth + 1 if node["type"] != "reaction" else depth
            dfs_traverse(child, new_depth)

    # Start traversal from the root
    dfs_traverse(route)

    print(
        f"Final result - alcohol_to_ketone: {alcohol_to_ketone}, ketone_to_alcohol: {ketone_to_alcohol}"
    )
    
    result = alcohol_to_ketone and ketone_to_alcohol

    # Add structural constraint if both transformations are found
    if result:
        structural_constraint_obj = {
            "type": "co-occurrence",
            "details": {
                "targets": [
                    "alcohol_to_carbonyl_oxidation",
                    "carbonyl_to_alcohol_reduction"
                ]
            }
        }
        if structural_constraint_obj not in findings_json["structural_constraints"]:
            findings_json["structural_constraints"].append(structural_constraint_obj)

    # Return True if both transformations are found
    return result, findings_json
