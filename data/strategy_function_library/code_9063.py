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
    Detects the use of ester formation as a protection strategy for carboxylic acids.
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
    result = False

    # Track protection and deprotection events
    protection_events = []
    deprotection_events = []

    def has_other_reactive_groups(smiles):
        """Check if molecule has other reactive functional groups that might need protection"""
        reactive_groups = [
            "Primary alcohol",
            "Secondary alcohol",
            "Tertiary alcohol",
            "Primary amine",
            "Secondary amine",
            "Aldehyde",
            "Ketone",
            "Thiol",
            "Phenol",
            "Alkyne",
            "Alkene",
        ]
        for group in reactive_groups:
            if checker.check_fg(group, smiles):
                if group not in findings_json["atomic_checks"]["functional_groups"]:
                    findings_json["atomic_checks"]["functional_groups"].append(group)
                return True
        return False

    def dfs_traverse(node, depth=0):
        nonlocal result
        if node["type"] == "reaction":
            if "metadata" in node and "mapped_reaction_smiles" in node["metadata"]:
                rsmi = node["metadata"]["mapped_reaction_smiles"]
                reactants = rsmi.split(">")[0].split(".")
                product = rsmi.split(">")[-1]

                # Check for esterification (protection)
                is_esterification = checker.check_reaction(
                    "Esterification of Carboxylic Acids", rsmi
                )

                if is_esterification:
                    if "Esterification of Carboxylic Acids" not in findings_json["atomic_checks"]["named_reactions"]:
                        findings_json["atomic_checks"]["named_reactions"].append("Esterification of Carboxylic Acids")
                    # Verify carboxylic acid in reactants and ester in product
                    for reactant in reactants:
                        if checker.check_fg("Carboxylic acid", reactant):
                            if "Carboxylic acid" not in findings_json["atomic_checks"]["functional_groups"]:
                                findings_json["atomic_checks"]["functional_groups"].append("Carboxylic acid")
                            if checker.check_fg("Ester", product):
                                if "Ester" not in findings_json["atomic_checks"]["functional_groups"]:
                                    findings_json["atomic_checks"]["functional_groups"].append("Ester")
                                protection_events.append((depth, reactant, product))

                # Check for deprotection
                is_deprotection = False
                deprotection_reactions = [
                    "Ester saponification (methyl deprotection)",
                    "Ester saponification (alkyl deprotection)",
                    "Hydrolysis or Hydrogenolysis of Carboxylic Esters or Thioesters"
                ]
                for rxn_name in deprotection_reactions:
                    if checker.check_reaction(rxn_name, rsmi):
                        is_deprotection = True
                        if rxn_name not in findings_json["atomic_checks"]["named_reactions"]:
                            findings_json["atomic_checks"]["named_reactions"].append(rxn_name)
                        break

                if is_deprotection:
                    # Verify ester in reactants and carboxylic acid in product
                    for reactant in reactants:
                        if checker.check_fg("Ester", reactant):
                            if "Ester" not in findings_json["atomic_checks"]["functional_groups"]:
                                findings_json["atomic_checks"]["functional_groups"].append("Ester")
                            if checker.check_fg("Carboxylic acid", product):
                                if "Carboxylic acid" not in findings_json["atomic_checks"]["functional_groups"]:
                                    findings_json["atomic_checks"]["functional_groups"].append("Carboxylic acid")
                                deprotection_events.append((depth, reactant, product))

        for child in node.get("children", []):
            # New logic for depth calculation
            new_depth = depth
            if node["type"] != "reaction": # If current node is chemical, depth increases
                new_depth = depth + 1
            # If current node is reaction, depth remains the same
            dfs_traverse(child, new_depth)

    dfs_traverse(route)

    # Analyze if protection strategy was used
    # Case 1: Both protection and deprotection events exist
    if protection_events and deprotection_events:
        # In retrosynthesis, protection should have higher depth than deprotection
        # (protection happens earlier in forward synthesis)
        for prot_depth, prot_reactant, prot_product in protection_events:
            for deprot_depth, deprot_reactant, deprot_product in deprotection_events:
                if prot_depth > deprot_depth:
                    result = True
                    # Add structural constraint for sequence
                    constraint = {
                        "type": "sequence",
                        "details": {
                            "before": "Esterification of Carboxylic Acids",
                            "after": [
                                "Ester saponification (methyl deprotection)",
                                "Ester saponification (alkyl deprotection)",
                                "Hydrolysis or Hydrogenolysis of Carboxylic Esters or Thioesters"
                            ]
                        }
                    }
                    if constraint not in findings_json["structural_constraints"]:
                        findings_json["structural_constraints"].append(constraint)
                    # No need to break, continue to find all instances

    # Case 2: Only protection events exist
    if protection_events:
        # Check if these protection events are likely for protection rather than the main goal
        for depth, reactant, product in protection_events:
            other_reactive_groups_in_reactant = has_other_reactive_groups(reactant)
            other_reactive_groups_in_product = has_other_reactive_groups(product)

            if other_reactive_groups_in_reactant or other_reactive_groups_in_product:
                result = True
                # Add structural constraints for co-occurrence with other reactive groups
                reactive_groups_list = [
                    "Primary alcohol", "Secondary alcohol", "Tertiary alcohol",
                    "Primary amine", "Secondary amine", "Aldehyde", "Ketone",
                    "Thiol", "Phenol", "Alkyne", "Alkene"
                ]
                for group in reactive_groups_list:
                    if checker.check_fg(group, reactant) or checker.check_fg(group, product):
                        constraint = {
                            "type": "co-occurrence",
                            "details": {
                                "targets": [
                                    "Esterification of Carboxylic Acids",
                                    group
                                ]
                            }
                        }
                        if constraint not in findings_json["structural_constraints"]:
                            findings_json["structural_constraints"].append(constraint)

    return result, findings_json
