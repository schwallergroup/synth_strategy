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


AMIDE_FORMATION_REACTIONS = [
    "Acylation of Nitrogen Nucleophiles by Acyl/Thioacyl/Carbamoyl Halides and Analogs_N",
    "Acylation of Nitrogen Nucleophiles by Carboxylic Acids",
    "Acyl chloride with secondary amine to amide",
    "Carboxylic acid with primary amine to amide",
    "Ester with secondary amine to amide",
]

def main(route) -> Tuple[bool, Dict]:
    """
    Detects a specific reaction sequence for nitrogen elaboration:
    (aldehyde or ketone) → oxime → primary amine → secondary amine → amide
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

    # Track the depth at which each transformation is found
    transformation_depths = {
        "oxime_formation": -1,
        "oxime_reduction": -1,
        "reductive_amination": -1,
        "amide_formation": -1,
    }

    def dfs_traverse(node, depth=0):
        nonlocal findings_json
        if node["type"] == "reaction":
            if "metadata" in node and "mapped_reaction_smiles" in node["metadata"]:
                rsmi = node["metadata"]["mapped_reaction_smiles"]
                reactants = rsmi.split(">")[0].split(".")
                product = rsmi.split(">")[-1]

                # Check for oxime formation (aldehyde/ketone + hydroxylamine → oxime)
                aldehyde_found = any(checker.check_fg("Aldehyde", r) for r in reactants)
                ketone_found = any(checker.check_fg("Ketone", r) for r in reactants)
                oxime_product_found = checker.check_fg("Oxime", product)

                if aldehyde_found or ketone_found:
                    if aldehyde_found and "Aldehyde" not in findings_json["atomic_checks"]["functional_groups"]:
                        findings_json["atomic_checks"]["functional_groups"].append("Aldehyde")
                    if ketone_found and "Ketone" not in findings_json["atomic_checks"]["functional_groups"]:
                        findings_json["atomic_checks"]["functional_groups"].append("Ketone")

                if oxime_product_found:
                    if "Oxime" not in findings_json["atomic_checks"]["functional_groups"]:
                        findings_json["atomic_checks"]["functional_groups"].append("Oxime")

                if (aldehyde_found or ketone_found) and oxime_product_found:
                    print(f"Found oxime formation at depth {depth}, rsmi: {rsmi}")
                    if transformation_depths["oxime_formation"] == -1:
                        transformation_depths["oxime_formation"] = depth
                    if "oxime_formation" not in findings_json["atomic_checks"]["named_reactions"]:
                        findings_json["atomic_checks"]["named_reactions"].append("oxime_formation")

                # Check for oxime reduction (oxime → primary amine)
                oxime_reactant_found = any(checker.check_fg("Oxime", r) for r in reactants)
                primary_amine_product_found = checker.check_fg("Primary amine", product)

                if oxime_reactant_found:
                    if "Oxime" not in findings_json["atomic_checks"]["functional_groups"]:
                        findings_json["atomic_checks"]["functional_groups"].append("Oxime")
                if primary_amine_product_found:
                    if "Primary amine" not in findings_json["atomic_checks"]["functional_groups"]:
                        findings_json["atomic_checks"]["functional_groups"].append("Primary amine")

                if oxime_reactant_found and primary_amine_product_found:
                    print(f"Found oxime reduction at depth {depth}, rsmi: {rsmi}")
                    if transformation_depths["oxime_reduction"] == -1:
                        transformation_depths["oxime_reduction"] = depth
                    if "oxime_reduction" not in findings_json["atomic_checks"]["named_reactions"]:
                        findings_json["atomic_checks"]["named_reactions"].append("oxime_reduction")

                # Check for reductive amination (aldehyde/ketone + primary amine → secondary amine)
                primary_amine_reactant_found = any(checker.check_fg("Primary amine", r) for r in reactants)
                secondary_amine_product_found = checker.check_fg("Secondary amine", product)

                if aldehyde_found or ketone_found:
                    if aldehyde_found and "Aldehyde" not in findings_json["atomic_checks"]["functional_groups"]:
                        findings_json["atomic_checks"]["functional_groups"].append("Aldehyde")
                    if ketone_found and "Ketone" not in findings_json["atomic_checks"]["functional_groups"]:
                        findings_json["atomic_checks"]["functional_groups"].append("Ketone")
                if primary_amine_reactant_found:
                    if "Primary amine" not in findings_json["atomic_checks"]["functional_groups"]:
                        findings_json["atomic_checks"]["functional_groups"].append("Primary amine")
                if secondary_amine_product_found:
                    if "Secondary amine" not in findings_json["atomic_checks"]["functional_groups"]:
                        findings_json["atomic_checks"]["functional_groups"].append("Secondary amine")

                if (
                    (aldehyde_found or ketone_found)
                    and primary_amine_reactant_found
                    and secondary_amine_product_found
                ):
                    print(f"Found reductive amination at depth {depth}, rsmi: {rsmi}")
                    if transformation_depths["reductive_amination"] == -1:
                        transformation_depths["reductive_amination"] = depth
                    if "reductive_amination" not in findings_json["atomic_checks"]["named_reactions"]:
                        findings_json["atomic_checks"]["named_reactions"].append("reductive_amination")

                # Check for amide formation using a list of known reaction types
                for reaction_name in AMIDE_FORMATION_REACTIONS:
                    if checker.check_reaction(reaction_name, rsmi):
                        print(f"Found amide formation at depth {depth}, rsmi: {rsmi}")
                        if transformation_depths["amide_formation"] == -1:
                            transformation_depths["amide_formation"] = depth
                        if reaction_name not in findings_json["atomic_checks"]["named_reactions"]:
                            findings_json["atomic_checks"]["named_reactions"].append(reaction_name)
                        break # Only need to find one amide formation reaction

        # Recursively process children with updated depth logic
        for child in node.get("children", []):
            new_depth = depth
            if node["type"] != "reaction": # Only increase depth when going from chemical to reaction
                new_depth = depth + 1
            dfs_traverse(child, new_depth)

    # Start traversal
    dfs_traverse(route)

    # Check if all transformations were found
    all_found = all(depth != -1 for depth in transformation_depths.values())

    # Check if they occurred in the correct sequence
    # Note: In retrosynthesis, higher depth = earlier in the synthesis
    correct_sequence = False
    if all_found:
        # The sequence should be: aldehyde → oxime → primary amine → secondary amine → amide
        # In retrosynthesis, this means: amide formation (latest) → reductive amination → oxime reduction → oxime formation (earliest)
        correct_sequence = (
            transformation_depths["amide_formation"]
            < transformation_depths["reductive_amination"]
            < transformation_depths["oxime_reduction"]
            < transformation_depths["oxime_formation"]
        )

    print(f"Transformation depths: {transformation_depths}")
    print(f"All transformations found: {all_found}")
    print(f"Correct sequence: {correct_sequence}")

    result = all_found and correct_sequence

    # Add structural constraints if conditions are met
    if all_found:
        # Co-occurrence constraint
        findings_json["structural_constraints"].append({
            "type": "co-occurrence",
            "details": {
                "targets": [
                    "oxime_formation",
                    "oxime_reduction",
                    "reductive_amination",
                    "amide_formation"
                ],
                "description": "All four key transformations (oxime formation, oxime reduction, reductive amination, and amide formation) must be present in the synthesis route."
            }
        })
    if correct_sequence:
        # Sequence constraint
        findings_json["structural_constraints"].append({
            "type": "sequence",
            "details": {
                "ordered_events": [
                    "oxime_formation",
                    "oxime_reduction",
                    "reductive_amination",
                    "amide_formation"
                ],
                "description": "The transformations must occur in a specific order, starting with oxime formation and ending with amide formation."
            }
        })

    return result, findings_json
