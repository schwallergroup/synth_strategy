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
    Detects the use of an ester as a protecting group for a carboxylic acid.
    This strategy is identified by finding an acid-to-ester transformation (protection)
    at an early synthetic stage, followed by an ester-to-acid transformation (deprotection)
    at a later stage.
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

    # Track the sequence of functional group transformations
    transformations = []

    def dfs_traverse(node, depth=0):
        nonlocal transformations, findings_json
        if node["type"] == "reaction":
            if "mapped_reaction_smiles" in node["metadata"]:
                rsmi = node["metadata"]["mapped_reaction_smiles"]
                reactants = rsmi.split(">")[0].split(".")
                product = rsmi.split(">")[-1]

                # Direct check for functional groups
                has_ester_reactant = any(checker.check_fg("Ester", r) for r in reactants)
                if has_ester_reactant:
                    if "Ester" not in findings_json["atomic_checks"]["functional_groups"]:
                        findings_json["atomic_checks"]["functional_groups"].append("Ester")

                has_acid_reactant = any(checker.check_fg("Carboxylic acid", r) for r in reactants)
                if has_acid_reactant:
                    if "Carboxylic acid" not in findings_json["atomic_checks"]["functional_groups"]:
                        findings_json["atomic_checks"]["functional_groups"].append("Carboxylic acid")

                has_ester_product = checker.check_fg("Ester", product)
                if has_ester_product:
                    if "Ester" not in findings_json["atomic_checks"]["functional_groups"]:
                        findings_json["atomic_checks"]["functional_groups"].append("Ester")

                has_acid_product = checker.check_fg("Carboxylic acid", product)
                if has_acid_product:
                    if "Carboxylic acid" not in findings_json["atomic_checks"]["functional_groups"]:
                        findings_json["atomic_checks"]["functional_groups"].append("Carboxylic acid")

                # Check for ester to acid transformation (hydrolysis/deprotection)
                deprotection_reactions = [
                    "Hydrolysis or Hydrogenolysis of Carboxylic Esters or Thioesters",
                    "Ester saponification (methyl deprotection)",
                    "Ester saponification (alkyl deprotection)"
                ]
                deprotection_found = False
                for r_name in deprotection_reactions:
                    if checker.check_reaction(r_name, rsmi):
                        if r_name not in findings_json["atomic_checks"]["named_reactions"]:
                            findings_json["atomic_checks"]["named_reactions"].append(r_name)
                        deprotection_found = True

                if deprotection_found:
                    # Verify ester in reactants and acid in product
                    if has_ester_reactant and has_acid_product:
                        transformations.append(("ester_to_acid", depth))
                        # Record structural constraint for deprotection event
                        if {"type": "co-occurrence", "details": {"targets": ["Hydrolysis or Hydrogenolysis of Carboxylic Esters or Thioesters", "Ester saponification (methyl deprotection)", "Ester saponification (alkyl deprotection)", "Ester", "Carboxylic acid"], "description": "Defines the 'ester_to_acid' deprotection event. One of the listed reactions must occur, consuming an 'Ester' reactant and producing a 'Carboxylic acid' product."}} not in findings_json["structural_constraints"]:
                            findings_json["structural_constraints"].append({"type": "co-occurrence", "details": {"targets": ["Hydrolysis or Hydrogenolysis of Carboxylic Esters or Thioesters", "Ester saponification (methyl deprotection)", "Ester saponification (alkyl deprotection)", "Ester", "Carboxylic acid"], "description": "Defines the 'ester_to_acid' deprotection event. One of the listed reactions must occur, consuming an 'Ester' reactant and producing a 'Carboxylic acid' product."}})

                # Check for acid to ester transformation (esterification/protection)
                protection_reactions = [
                    "Esterification of Carboxylic Acids",
                    "O-alkylation of carboxylic acids with diazo compounds",
                    "Transesterification"
                ]
                protection_found = False
                for r_name in protection_reactions:
                    if checker.check_reaction(r_name, rsmi):
                        if r_name not in findings_json["atomic_checks"]["named_reactions"]:
                            findings_json["atomic_checks"]["named_reactions"].append(r_name)
                        protection_found = True

                if protection_found:
                    # Verify acid in reactants and ester in product
                    if has_acid_reactant and has_ester_product:
                        transformations.append(("acid_to_ester", depth))
                        # Record structural constraint for protection event
                        if {"type": "co-occurrence", "details": {"targets": ["Esterification of Carboxylic Acids", "O-alkylation of carboxylic acids with diazo compounds", "Transesterification", "Carboxylic acid", "Ester"], "description": "Defines the 'acid_to_ester' protection event. One of the listed reactions must occur, consuming a 'Carboxylic acid' reactant and producing an 'Ester' product."}} not in findings_json["structural_constraints"]:
                            findings_json["structural_constraints"].append({"type": "co-occurrence", "details": {"targets": ["Esterification of Carboxylic Acids", "O-alkylation of carboxylic acids with diazo compounds", "Transesterification", "Carboxylic acid", "Ester"], "description": "Defines the 'acid_to_ester' protection event. One of the listed reactions must occur, consuming a 'Carboxylic acid' reactant and producing an 'Ester' product."}})

        # Traverse children
        for child in node.get("children", []):
            # New logic for depth calculation
            new_depth = depth
            if node["type"] != "reaction": # If current node is chemical, depth increases
                new_depth = depth + 1
            # If current node is reaction, depth remains the same
            dfs_traverse(child, new_depth)

    # Start traversal from the root
    dfs_traverse(route)

    # Check if we have the correct protection -> deprotection sequence
    # In synthesis, higher depth means earlier stage.
    has_sequence = False
    for i in range(len(transformations)):
        for j in range(len(transformations)):
            # Look for protection (acid->ester) at an earlier step (higher depth)
            # and deprotection (ester->acid) at a later step (lower depth).
            if (
                transformations[i][0] == "acid_to_ester"
                and transformations[j][0] == "ester_to_acid"
                and transformations[i][1] > transformations[j][1]
            ):
                has_sequence = True
                # Record structural constraint for sequence
                if {"type": "sequence", "details": {"ordered_events": ["acid_to_ester_transformation", "ester_to_acid_transformation"], "description": "The 'acid_to_ester' transformation must occur at an earlier synthetic stage (higher depth value) than the 'ester_to_acid' transformation."}} not in findings_json["structural_constraints"]:
                    findings_json["structural_constraints"].append({"type": "sequence", "details": {"ordered_events": ["acid_to_ester_transformation", "ester_to_acid_transformation"], "description": "The 'acid_to_ester' transformation must occur at an earlier synthetic stage (higher depth value) than the 'ester_to_acid' transformation."}})
                break
        if has_sequence:
            break

    return has_sequence, findings_json
