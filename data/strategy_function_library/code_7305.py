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

def main(route) -> Tuple[bool, Dict]:
    """
    Detects if the synthesis uses a late-stage Grignard addition to an acid chloride
    to form a C-C bond, following ester hydrolysis and acid chloride formation.
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

    # Track if we found each step in our sequence
    found_ester_hydrolysis = False
    found_acid_chloride_formation = False
    found_grignard_addition = False

    # Track depths for determining stage of synthesis
    max_depth = 0
    grignard_depth = None

    # Track reaction sequence
    reaction_sequence = []

    def dfs_traverse(node, depth=0):
        nonlocal found_ester_hydrolysis, found_acid_chloride_formation, found_grignard_addition
        nonlocal max_depth, grignard_depth, reaction_sequence, findings_json

        max_depth = max(max_depth, depth)

        if node["type"] == "reaction" and "metadata" in node and "mapped_reaction_smiles" in node["metadata"]:
            rsmi = node["metadata"]["mapped_reaction_smiles"]
            reactants = rsmi.split(">")[0].split(".")
            product = rsmi.split(">")[-1]

            # Check for ester hydrolysis: ester to carboxylic acid
            if (
                checker.check_reaction("Ester saponification (alkyl deprotection)", rsmi)
                or checker.check_reaction("Ester saponification (methyl deprotection)", rsmi)
            ):
                found_ester_hydrolysis = True
                reaction_sequence.append(("ester_hydrolysis", depth))
                if "Ester saponification (alkyl deprotection)" not in findings_json["atomic_checks"]["named_reactions"]:
                    findings_json["atomic_checks"]["named_reactions"].append("Ester saponification (alkyl deprotection)")
                if "Ester saponification (methyl deprotection)" not in findings_json["atomic_checks"]["named_reactions"]:
                    findings_json["atomic_checks"]["named_reactions"].append("Ester saponification (methyl deprotection)")
                print(f"Found ester hydrolysis at depth {depth}")

            # Check for acid chloride formation: carboxylic acid to acid chloride
            if (
                any(checker.check_fg("Carboxylic acid", r) for r in reactants)
                and checker.check_fg("Acyl halide", product)
            ):
                found_acid_chloride_formation = True
                reaction_sequence.append(("acid_chloride_formation", depth))
                if "Carboxylic acid" not in findings_json["atomic_checks"]["functional_groups"]:
                    findings_json["atomic_checks"]["functional_groups"].append("Carboxylic acid")
                if "Acyl halide" not in findings_json["atomic_checks"]["functional_groups"]:
                    findings_json["atomic_checks"]["functional_groups"].append("Acyl halide")
                print(f"Found acid chloride formation at depth {depth}")

            # Check for Grignard addition: acid chloride + Grignard reagent to ketone
            if (
                any(checker.check_fg("Magnesium halide", r) for r in reactants)
                and any(checker.check_fg("Acyl halide", r) for r in reactants)
                and checker.check_fg("Ketone", product)
            ):
                found_grignard_addition = True
                grignard_depth = depth
                reaction_sequence.append(("grignard_addition", depth))
                if "Magnesium halide" not in findings_json["atomic_checks"]["functional_groups"]:
                    findings_json["atomic_checks"]["functional_groups"].append("Magnesium halide")
                if "Acyl halide" not in findings_json["atomic_checks"]["functional_groups"]:
                    findings_json["atomic_checks"]["functional_groups"].append("Acyl halide")
                if "Ketone" not in findings_json["atomic_checks"]["functional_groups"]:
                    findings_json["atomic_checks"]["functional_groups"].append("Ketone")
                print(f"Found Grignard addition at depth {depth}")

        # Recursively process children
        for child in node.get("children", []):
            # New logic for depth calculation
            if node["type"] == "reaction":
                # If current node is reaction, depth remains same for child (chemical)
                dfs_traverse(child, depth)
            else:
                # If current node is chemical, depth increases for child (reaction)
                dfs_traverse(child, depth + 1)

    # Start traversal from the root
    dfs_traverse(route)

    # Check if we found all steps
    all_steps_found = (
        found_ester_hydrolysis and found_acid_chloride_formation and found_grignard_addition
    )

    # Check if Grignard addition is late-stage (in first 40% of synthesis depth)
    is_late_stage = grignard_depth is not None and grignard_depth <= max_depth * 0.4

    # Check if reactions are in the correct sequence for retrosynthetic traversal
    # In retrosynthesis, we expect: Grignard -> acid chloride -> ester hydrolysis
    correct_sequence = False
    if all_steps_found:
        # Sort reactions by depth (ascending)
        reaction_sequence.sort(key=lambda x: x[1])
        # Extract just the reaction types
        sequence_types = [r[0] for r in reaction_sequence]

        print(f"Sequence types in order: {sequence_types}")

        # Check if the sequence contains our pattern in retrosynthetic order
        for i in range(len(sequence_types) - 2):
            if (
                sequence_types[i] == "grignard_addition"
                and sequence_types[i + 1] == "acid_chloride_formation"
                and sequence_types[i + 2] == "ester_hydrolysis"
            ):
                correct_sequence = True
                break

    print(f"All steps found: {all_steps_found}")
    print(f"Is late stage: {is_late_stage}")
    print(f"Correct sequence: {correct_sequence}")
    print(f"Reaction sequence: {reaction_sequence}")

    result = all_steps_found and is_late_stage and correct_sequence

    # Populate structural constraints based on the final checks
    if all_steps_found:
        findings_json["structural_constraints"].append({
            "type": "co-occurrence",
            "details": {
                "targets": [
                    "ester_hydrolysis",
                    "acid_chloride_formation",
                    "grignard_addition"
                ],
                "description": "The route must contain an ester hydrolysis step, an acid chloride formation step, and a Grignard addition step."
            }
        })
    
    if is_late_stage:
        findings_json["structural_constraints"].append({
            "type": "positional",
            "details": {
                "target": "grignard_addition",
                "position": "depth <= 0.4 * max_depth",
                "description": "The Grignard addition must occur in the first 40% of retrosynthetic steps (i.e., late-stage in the forward synthesis)."
            }
        })

    if correct_sequence:
        findings_json["structural_constraints"].append({
            "type": "sequence",
            "details": {
                "ordered_events": [
                    "grignard_addition",
                    "acid_chloride_formation",
                    "ester_hydrolysis"
                ],
                "basis": "retrosynthetic_depth",
                "description": "The events must occur in a specific order based on their depth in the retrosynthetic tree."
            }
        })

    return result, findings_json
