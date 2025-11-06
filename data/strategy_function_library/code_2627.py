from typing import Tuple, Dict, List
import copy
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
    Detects if the synthesis route employs a protection-deprotection sequence,
    specifically looking for alcohol, amine, or carboxylic acid protection and deprotection.
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

    # Track protection and deprotection steps with their depths and functional groups
    protection_steps = []
    deprotection_steps = []

    def dfs_traverse(node, depth=0):
        nonlocal protection_steps, deprotection_steps, findings_json
        if node["type"] == "reaction" and "mapped_reaction_smiles" in node["metadata"]:
            rsmi = node["metadata"]["mapped_reaction_smiles"]
            reactants = rsmi.split(">")[0].split(".")
            product = rsmi.split(">")[-1]

            # Check for protection reactions using reaction checkers
            protection_found = False
            if checker.check_reaction("Alcohol protection with silyl ethers", rsmi):
                findings_json["atomic_checks"]["named_reactions"].append("Alcohol protection with silyl ethers")
                protection_found = True
            if checker.check_reaction("Boc amine protection", rsmi):
                findings_json["atomic_checks"]["named_reactions"].append("Boc amine protection")
                protection_found = True
            if checker.check_reaction("Protection of carboxylic acid", rsmi):
                findings_json["atomic_checks"]["named_reactions"].append("Protection of carboxylic acid")
                protection_found = True

            if protection_found:
                # Identify which functional group is being protected
                fg_type = None
                if any(
                    checker.check_fg("Primary alcohol", r)
                    or checker.check_fg("Secondary alcohol", r)
                    or checker.check_fg("Tertiary alcohol", r)
                    or checker.check_fg("Phenol", r)
                    for r in reactants
                ):
                    fg_type = "alcohol"
                    for r in reactants:
                        if checker.check_fg("Primary alcohol", r):
                            findings_json["atomic_checks"]["functional_groups"].append("Primary alcohol")
                        if checker.check_fg("Secondary alcohol", r):
                            findings_json["atomic_checks"]["functional_groups"].append("Secondary alcohol")
                        if checker.check_fg("Tertiary alcohol", r):
                            findings_json["atomic_checks"]["functional_groups"].append("Tertiary alcohol")
                        if checker.check_fg("Phenol", r):
                            findings_json["atomic_checks"]["functional_groups"].append("Phenol")
                elif any(
                    checker.check_fg("Primary amine", r) or checker.check_fg("Secondary amine", r)
                    for r in reactants
                ):
                    fg_type = "amine"
                    for r in reactants:
                        if checker.check_fg("Primary amine", r):
                            findings_json["atomic_checks"]["functional_groups"].append("Primary amine")
                        if checker.check_fg("Secondary amine", r):
                            findings_json["atomic_checks"]["functional_groups"].append("Secondary amine")
                elif any(checker.check_fg("Carboxylic acid", r) for r in reactants):
                    fg_type = "carboxylic_acid"
                    for r in reactants:
                        if checker.check_fg("Carboxylic acid", r):
                            findings_json["atomic_checks"]["functional_groups"].append("Carboxylic acid")

                if fg_type:
                    protection_steps.append((depth, fg_type))

            # Check for deprotection reactions using reaction checkers
            deprotection_found = False
            if checker.check_reaction("Alcohol deprotection from silyl ethers", rsmi):
                findings_json["atomic_checks"]["named_reactions"].append("Alcohol deprotection from silyl ethers")
                deprotection_found = True
            if checker.check_reaction("Boc amine deprotection", rsmi):
                findings_json["atomic_checks"]["named_reactions"].append("Boc amine deprotection")
                deprotection_found = True
            if checker.check_reaction("Deprotection of carboxylic acid", rsmi):
                findings_json["atomic_checks"]["named_reactions"].append("Deprotection of carboxylic acid")
                deprotection_found = True
            if checker.check_reaction("Hydrolysis or Hydrogenolysis of Carboxylic Esters or Thioesters", rsmi):
                findings_json["atomic_checks"]["named_reactions"].append("Hydrolysis or Hydrogenolysis of Carboxylic Esters or Thioesters")
                deprotection_found = True
            if checker.check_reaction("COOH ethyl deprotection", rsmi):
                findings_json["atomic_checks"]["named_reactions"].append("COOH ethyl deprotection")
                deprotection_found = True
            if checker.check_reaction("Ester saponification (methyl deprotection)", rsmi):
                findings_json["atomic_checks"]["named_reactions"].append("Ester saponification (methyl deprotection)")
                deprotection_found = True
            if checker.check_reaction("Ester saponification (alkyl deprotection)", rsmi):
                findings_json["atomic_checks"]["named_reactions"].append("Ester saponification (alkyl deprotection)")
                deprotection_found = True

            if deprotection_found:
                # Identify which functional group is being deprotected
                fg_type = None
                if (
                    checker.check_fg("Primary alcohol", product)
                    or checker.check_fg("Secondary alcohol", product)
                    or checker.check_fg("Tertiary alcohol", product)
                    or checker.check_fg("Phenol", product)
                ):
                    fg_type = "alcohol"
                    if checker.check_fg("Primary alcohol", product):
                        findings_json["atomic_checks"]["functional_groups"].append("Primary alcohol")
                    if checker.check_fg("Secondary alcohol", product):
                        findings_json["atomic_checks"]["functional_groups"].append("Secondary alcohol")
                    if checker.check_fg("Tertiary alcohol", product):
                        findings_json["atomic_checks"]["functional_groups"].append("Tertiary alcohol")
                    if checker.check_fg("Phenol", product):
                        findings_json["atomic_checks"]["functional_groups"].append("Phenol")
                elif checker.check_fg("Primary amine", product) or checker.check_fg(
                    "Secondary amine", product
                ):
                    fg_type = "amine"
                    if checker.check_fg("Primary amine", product):
                        findings_json["atomic_checks"]["functional_groups"].append("Primary amine")
                    if checker.check_fg("Secondary amine", product):
                        findings_json["atomic_checks"]["functional_groups"].append("Secondary amine")
                elif checker.check_fg("Carboxylic acid", product):
                    fg_type = "carboxylic_acid"
                    if checker.check_fg("Carboxylic acid", product):
                        findings_json["atomic_checks"]["functional_groups"].append("Carboxylic acid")

                if fg_type:
                    deprotection_steps.append((depth, fg_type))

        for child in node.get("children", []):
            # New logic for depth calculation
            new_depth = depth
            if node["type"] != "reaction":  # If current node is chemical, depth increases
                new_depth = depth + 1
            dfs_traverse(child, new_depth)

    dfs_traverse(route)

    result = False
    # Check if we have both protection and deprotection steps
    if protection_steps and deprotection_steps:
        # Check for matching protection-deprotection pairs
        for prot_depth, prot_fg in protection_steps:
            for deprot_depth, deprot_fg in deprotection_steps:
                # In retrosynthesis, higher depth = earlier stage in forward synthesis
                # So for a proper protection-deprotection sequence, the protection depth
                # should be greater than the deprotection depth
                if prot_depth > deprot_depth and prot_fg == deprot_fg:
                    result = True
                    # Add the structural constraint if a valid sequence is found
                    findings_json["structural_constraints"].append({
                        "type": "sequence",
                        "details": {
                            "before": {
                                "event_name": "protection",
                                "triggered_by_reactions": [
                                    "Alcohol protection with silyl ethers",
                                    "Boc amine protection",
                                    "Protection of carboxylic acid"
                                ]
                            },
                            "after": {
                                "event_name": "deprotection",
                                "triggered_by_reactions": [
                                    "Alcohol deprotection from silyl ethers",
                                    "Boc amine deprotection",
                                    "Deprotection of carboxylic acid",
                                    "Hydrolysis or Hydrogenolysis of Carboxylic Esters or Thioesters",
                                    "COOH ethyl deprotection",
                                    "Ester saponification (methyl deprotection)",
                                    "Ester saponification (alkyl deprotection)"
                                ]
                            },
                            "condition": "The functional group type (alcohol, amine, or carboxylic acid) associated with the 'before' event must be the same as the functional group type associated with the 'after' event."
                        }
                    })
                    # Break after finding the first valid sequence to avoid duplicate structural constraints
                    # if multiple pairs satisfy the condition, but the overall result is still True.
                    break 
            if result: # If a valid sequence is found, no need to check further protection steps
                break

    # Ensure unique entries in atomic_checks lists
    findings_json["atomic_checks"]["named_reactions"] = list(set(findings_json["atomic_checks"]["named_reactions"]))
    findings_json["atomic_checks"]["functional_groups"] = list(set(findings_json["atomic_checks"]["functional_groups"]))

    return result, findings_json
