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


from rdkit import Chem

# Assuming 'checker' is a pre-defined object with the necessary methods.

HETEROARYL_COUPLING_REACTIONS = [
    "Suzuki coupling with boronic acids",
    "Suzuki coupling with boronic esters",
    "Negishi coupling",
    "Stille reaction_aryl",
    "Heck terminal vinyl",
]

def main(route) -> Tuple[bool, Dict]:
    """
    Detects a multi-step protect-modify-deprotect-couple sequence. The function specifically identifies the following ordered pattern: 1. Ester formation from a carboxylic acid (Protection). 2. A heteroaryl coupling reaction from a defined list including Suzuki, Negishi, Stille, and Heck (Modification). 3. Ester hydrolysis to a carboxylic acid (Deprotection). 4. Amide formation from the carboxylic acid (Coupling).
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

    # Track reaction sequence
    reaction_sequence = []

    def identify_reaction_type(reactants, product, rxn_smiles):
        """Identify the type of reaction"""
        # Check for ester formation (protection)
        if any(
            checker.check_fg("Carboxylic acid", r) for r in reactants if Chem.MolFromSmiles(r)
        ) and checker.check_fg("Ester", product):
            print(f"Found protection reaction: {rxn_smiles}")
            findings_json["atomic_checks"]["named_reactions"].append("ester_formation")
            findings_json["atomic_checks"]["functional_groups"].append("Carboxylic acid")
            findings_json["atomic_checks"]["functional_groups"].append("Ester")
            return "protection"

        # Check for coupling reactions (modification)
        for r_name in HETEROARYL_COUPLING_REACTIONS:
            if checker.check_reaction(r_name, rxn_smiles):
                print(f"Found modification reaction: {rxn_smiles}")
                findings_json["atomic_checks"]["named_reactions"].append(r_name)
                return "modification"

        # Check for ester hydrolysis (deprotection)
        if any(
            checker.check_fg("Ester", r) for r in reactants if Chem.MolFromSmiles(r)
        ) and checker.check_fg("Carboxylic acid", product):
            print(f"Found deprotection reaction: {rxn_smiles}")
            findings_json["atomic_checks"]["named_reactions"].append("ester_hydrolysis")
            findings_json["atomic_checks"]["functional_groups"].append("Ester")
            findings_json["atomic_checks"]["functional_groups"].append("Carboxylic acid")
            return "deprotection"

        # Check for amide formation (coupling)
        if any(checker.check_fg("Carboxylic acid", r) for r in reactants if Chem.MolFromSmiles(r)):
            if (
                checker.check_fg("Primary amide", product)
                or checker.check_fg("Secondary amide", product)
                or checker.check_fg("Tertiary amide", product)
            ):
                print(f"Found coupling reaction: {rxn_smiles}")
                findings_json["atomic_checks"]["named_reactions"].append("amide_formation")
                findings_json["atomic_checks"]["functional_groups"].append("Carboxylic acid")
                if checker.check_fg("Primary amide", product):
                    findings_json["atomic_checks"]["functional_groups"].append("Primary amide")
                if checker.check_fg("Secondary amide", product):
                    findings_json["atomic_checks"]["functional_groups"].append("Secondary amide")
                if checker.check_fg("Tertiary amide", product):
                    findings_json["atomic_checks"]["functional_groups"].append("Tertiary amide")
                return "coupling"

        print(f"Unidentified reaction type: {rxn_smiles}")
        return "other"

    def dfs_traverse(node, depth=0):
        if node["type"] == "reaction":
            if "metadata" in node and "mapped_reaction_smiles" in node["metadata"]:
                rsmi = node["metadata"]["mapped_reaction_smiles"]
                reactants = rsmi.split(">")[0].split(".")
                product = rsmi.split(">")[-1]

                reaction_type = identify_reaction_type(reactants, product, rsmi)
                reaction_sequence.append((reaction_type, depth, rsmi))
                print(f"Depth {depth}: Identified reaction: {reaction_type}")

        for child in node.get("children", []):
            # New logic for depth calculation
            new_depth = depth
            if node["type"] != "reaction": # This means current node is 'chemical'
                new_depth = depth + 1
            dfs_traverse(child, new_depth)

    # Start traversal
    dfs_traverse(route)

    # Sort by depth to get synthetic (forward) direction
    reaction_sequence.sort(key=lambda x: x[1], reverse=True)
    forward_sequence = [r[0] for r in reaction_sequence]

    print(f"Reaction sequence (forward direction): {forward_sequence}")

    # Check if the pattern exists in the sequence
    has_protection = "protection" in forward_sequence
    has_modification = "modification" in forward_sequence
    has_deprotection = "deprotection" in forward_sequence
    has_coupling = "coupling" in forward_sequence

    # Check if they appear in the correct order
    correct_order = False
    if has_protection and has_modification and has_deprotection and has_coupling:
        protection_idx = forward_sequence.index("protection")
        modification_indices = [i for i, x in enumerate(forward_sequence) if x == "modification"]
        deprotection_idx = forward_sequence.index("deprotection")
        coupling_idx = forward_sequence.index("coupling")

        # Check if at least one modification is between protection and deprotection
        modification_between = any(
            protection_idx < mod_idx < deprotection_idx for mod_idx in modification_indices
        )

        # Check if deprotection is before coupling
        deprotection_before_coupling = deprotection_idx < coupling_idx

        correct_order = modification_between and deprotection_before_coupling

    strategy_present = (
        has_protection and has_modification and has_deprotection and has_coupling and correct_order
    )

    print(f"Strategy detection results:")
    print(f"- Has protection: {has_protection}")
    print(f"- Has modification: {has_modification}")
    print(f"- Has deprotection: {has_deprotection}")
    print(f"- Has coupling: {has_coupling}")
    print(f"- Correct order: {correct_order}")
    print(f"- Overall strategy present: {strategy_present}")

    if has_protection and has_modification and has_deprotection and has_coupling:
        findings_json["structural_constraints"].append({
            "type": "co-occurrence",
            "details": {
                "targets": [
                    "protection",
                    "modification",
                    "deprotection",
                    "coupling"
                ]
            }
        })
    if correct_order:
        findings_json["structural_constraints"].append({
            "type": "sequence",
            "details": {
                "ordered_steps": [
                    "protection",
                    "modification",
                    "deprotection",
                    "coupling"
                ]
            }
        })

    return strategy_present, findings_json
