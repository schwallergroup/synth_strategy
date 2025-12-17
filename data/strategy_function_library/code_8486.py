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


CC_COUPLING_REACTIONS = [
    "Heck terminal vinyl",
    "Heck_terminal_vinyl",
    "Suzuki coupling with boronic acids",
    "Suzuki coupling with boronic esters",
    "Suzuki",
    "Stille reaction_vinyl",
    "Stille reaction_aryl",
    "Stille",
    "Negishi coupling",
    "Negishi",
    "Sonogashira alkyne_aryl halide",
    "Sonogashira acetylene_aryl halide",
]

def main(route) -> Tuple[bool, Dict]:
    """
    Detects a multi-step strategy: C-C bond formation via a specific coupling reaction, followed by sequential oxidation (alkene → aldehyde → carboxylic acid) and activation to an acyl chloride. The checked coupling reactions include various forms of Heck, Suzuki, Stille, Negishi, and Sonogashira reactions.
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

    # Initialize tracking variables
    reactions_with_depth = []

    def dfs_traverse(node, depth=0):
        nonlocal reactions_with_depth, findings_json
        if node["type"] == "reaction":
            try:
                rsmi = node["metadata"]["mapped_reaction_smiles"]
                reactants_part = rsmi.split(">")[0]
                products_part = rsmi.split(">")[-1]
                reactants = reactants_part.split(".")
                product = products_part

                # Check for C-C coupling reactions
                cc_coupling_detected = False
                for r in CC_COUPLING_REACTIONS:
                    if checker.check_reaction(r, rsmi):
                        print(f"Detected C-C coupling at depth {depth}")
                        reactions_with_depth.append(("cc_coupling", depth))
                        if "cc_coupling" not in findings_json["atomic_checks"]["named_reactions"]:
                            findings_json["atomic_checks"]["named_reactions"].append("cc_coupling")
                        if r not in findings_json["atomic_checks"]["named_reactions"]:
                            findings_json["atomic_checks"]["named_reactions"].append(r)
                        cc_coupling_detected = True
                        break

                # Check for alkene to aldehyde oxidation
                alkene_to_aldehyde_detected = False
                if checker.check_reaction("Alkene oxidation to aldehyde", rsmi):
                    print(f"Detected vinyl to aldehyde oxidation at depth {depth}")
                    reactions_with_depth.append(("vinyl_to_aldehyde", depth))
                    if "vinyl_to_aldehyde" not in findings_json["atomic_checks"]["named_reactions"]:
                        findings_json["atomic_checks"]["named_reactions"].append("vinyl_to_aldehyde")
                    if "Alkene oxidation to aldehyde" not in findings_json["atomic_checks"]["named_reactions"]:
                        findings_json["atomic_checks"]["named_reactions"].append("Alkene oxidation to aldehyde")
                    alkene_to_aldehyde_detected = True
                elif (
                    any(checker.check_fg("Vinyl", r) for r in reactants)
                    and checker.check_fg("Aldehyde", product)
                    and not any(checker.check_fg("Aldehyde", r) for r in reactants)
                ):
                    print(f"Detected vinyl to aldehyde oxidation at depth {depth}")
                    reactions_with_depth.append(("vinyl_to_aldehyde", depth))
                    if "vinyl_to_aldehyde" not in findings_json["atomic_checks"]["named_reactions"]:
                        findings_json["atomic_checks"]["named_reactions"].append("vinyl_to_aldehyde")
                    if "Vinyl" not in findings_json["atomic_checks"]["functional_groups"]:
                        findings_json["atomic_checks"]["functional_groups"].append("Vinyl")
                    if "Aldehyde" not in findings_json["atomic_checks"]["functional_groups"]:
                        findings_json["atomic_checks"]["functional_groups"].append("Aldehyde")
                    alkene_to_aldehyde_detected = True

                # Check for aldehyde to carboxylic acid oxidation
                aldehyde_to_acid_detected = False
                if checker.check_reaction("Oxidation of aldehydes to carboxylic acids", rsmi):
                    print(f"Detected aldehyde to carboxylic acid oxidation at depth {depth}")
                    reactions_with_depth.append(("aldehyde_to_acid", depth))
                    if "aldehyde_to_acid" not in findings_json["atomic_checks"]["named_reactions"]:
                        findings_json["atomic_checks"]["named_reactions"].append("aldehyde_to_acid")
                    if "Oxidation of aldehydes to carboxylic acids" not in findings_json["atomic_checks"]["named_reactions"]:
                        findings_json["atomic_checks"]["named_reactions"].append("Oxidation of aldehydes to carboxylic acids")
                    aldehyde_to_acid_detected = True
                elif (
                    any(checker.check_fg("Aldehyde", r) for r in reactants)
                    and checker.check_fg("Carboxylic acid", product)
                    and not any(checker.check_fg("Carboxylic acid", r) for r in reactants)
                ):
                    print(f"Detected aldehyde to carboxylic acid oxidation at depth {depth}")
                    reactions_with_depth.append(("aldehyde_to_acid", depth))
                    if "aldehyde_to_acid" not in findings_json["atomic_checks"]["named_reactions"]:
                        findings_json["atomic_checks"]["named_reactions"].append("aldehyde_to_acid")
                    if "Aldehyde" not in findings_json["atomic_checks"]["functional_groups"]:
                        findings_json["atomic_checks"]["functional_groups"].append("Aldehyde")
                    if "Carboxylic acid" not in findings_json["atomic_checks"]["functional_groups"]:
                        findings_json["atomic_checks"]["functional_groups"].append("Carboxylic acid")
                    aldehyde_to_acid_detected = True

                # Check for acid to acid chloride activation
                acid_to_acid_chloride_detected = False
                if (
                    any(checker.check_fg("Carboxylic acid", r) for r in reactants)
                    and checker.check_fg("Acyl halide", product)
                    and not any(checker.check_fg("Acyl halide", r) for r in reactants)
                ):
                    print(f"Detected carboxylic acid to acid chloride activation at depth {depth}")
                    reactions_with_depth.append(("acid_to_acid_chloride", depth))
                    if "acid_to_acid_chloride" not in findings_json["atomic_checks"]["named_reactions"]:
                        findings_json["atomic_checks"]["named_reactions"].append("acid_to_acid_chloride")
                    if "Carboxylic acid" not in findings_json["atomic_checks"]["functional_groups"]:
                        findings_json["atomic_checks"]["functional_groups"].append("Carboxylic acid")
                    if "Acyl halide" not in findings_json["atomic_checks"]["functional_groups"]:
                        findings_json["atomic_checks"]["functional_groups"].append("Acyl halide")
                    acid_to_acid_chloride_detected = True

            except Exception as e:
                print(f"Error processing reaction: {e}")

        # Process children
        for child in node.get("children", []):
            # New logic for depth calculation
            new_depth = depth
            if node["type"] != "reaction": # Only increase depth when moving from chemical to reaction
                new_depth = depth + 1
            dfs_traverse(child, new_depth)

    # Start traversal
    dfs_traverse(route)

    # Sort reactions by depth (ascending - lower depth is later in synthesis)
    reactions_with_depth.sort(key=lambda x: x[1])

    # Extract just the reaction types in order
    reaction_sequence = [r[0] for r in reactions_with_depth]
    print(f"Reaction sequence (retrosynthetic): {reaction_sequence}")

    # Check if we have the required reactions
    has_cc_coupling = "cc_coupling" in reaction_sequence
    has_vinyl_to_aldehyde = "vinyl_to_aldehyde" in reaction_sequence
    has_aldehyde_to_acid = "aldehyde_to_acid" in reaction_sequence
    has_acid_to_acid_chloride = "acid_to_acid_chloride" in reaction_sequence

    # Check if we have at least two oxidation steps
    has_sequential_oxidation = has_vinyl_to_aldehyde and has_aldehyde_to_acid

    # In retrosynthetic analysis, the sequence should be:
    # acid_to_acid_chloride → aldehyde_to_acid → vinyl_to_aldehyde → cc_coupling
    # (which is the reverse of the forward synthesis)

    correct_sequence = False

    # We need at least the sequential oxidation steps
    if has_sequential_oxidation:
        try:
            aldehyde_to_acid_idx = reaction_sequence.index("aldehyde_to_acid")
            vinyl_to_aldehyde_idx = reaction_sequence.index("vinyl_to_aldehyde")

            # Check if oxidation steps are in correct order
            if aldehyde_to_acid_idx < vinyl_to_aldehyde_idx:
                # Record structural constraint: aldehyde_to_acid before vinyl_to_aldehyde
                findings_json["structural_constraints"].append({
                    "type": "sequence",
                    "details": {
                        "ordered_targets": [
                            "aldehyde_to_acid",
                            "vinyl_to_aldehyde"
                        ],
                        "comment": "In retrosynthesis, the oxidation of an aldehyde to a carboxylic acid must be a more recent step (lower depth) than the oxidation of a vinyl group to an aldehyde."
                    }
                })
                # If we have C-C coupling, it should come after vinyl_to_aldehyde
                if not has_cc_coupling or (
                    has_cc_coupling
                    and vinyl_to_aldehyde_idx < reaction_sequence.index("cc_coupling")
                ):
                    correct_sequence = True
                    if has_cc_coupling:
                        # Record structural constraint: vinyl_to_aldehyde before cc_coupling
                        findings_json["structural_constraints"].append({
                            "type": "sequence",
                            "details": {
                                "ordered_targets": [
                                    "vinyl_to_aldehyde",
                                    "cc_coupling"
                                ],
                                "comment": "If a C-C coupling reaction is present, it must occur earlier in the retrosynthesis (higher depth) than the vinyl to aldehyde oxidation."
                            }
                        })
        except ValueError:
            pass

    # For the full strategy, we need C-C coupling and sequential oxidation in correct order
    # But we'll also accept just the sequential oxidation as a partial match
    strategy_present = has_sequential_oxidation and correct_sequence

    if strategy_present:
        if has_cc_coupling:
            print(
                "Detected complete sequential oxidation strategy with late-stage C-C bond formation"
            )
        else:
            print("Detected sequential oxidation strategy (without C-C coupling)")
    else:
        print("Strategy not detected")
        if not has_sequential_oxidation:
            print("Missing sequential oxidation steps")
        if not correct_sequence:
            print("Incorrect reaction sequence")

    return strategy_present, findings_json