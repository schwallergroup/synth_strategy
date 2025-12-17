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


HETEROCYCLES_OF_INTEREST = [
    "oxazole",
    "thiazole",
    "imidazole",
    "benzoxazole",
    "benzothiazole",
    "benzimidazole",
]

def main(route) -> Tuple[bool, Dict]:
    """This function detects a linear synthesis strategy involving acid chloride formation, amide coupling, formation of a specific heterocycle, and late-stage SNAr diversification. The heterocycles of interest are defined in the HETEROCYCLES_OF_INTEREST list, including: oxazole, thiazole, imidazole, benzoxazole, benzothiazole, and benzimidazole."""
    findings_template = {
        "atomic_checks": {
            "named_reactions": [],
            "ring_systems": [],
            "functional_groups": []
        },
        "structural_constraints": []
    }
    findings_json = copy.deepcopy(findings_template)

    # Track each key step
    acid_chloride_reactions = []
    amide_formation_reactions = []
    heterocycle_formation_reactions = []
    snar_diversification_reactions = []

    # Store reaction sequence for verification
    reaction_sequence = []

    def dfs_traverse(node, depth=0):
        nonlocal acid_chloride_reactions, amide_formation_reactions, heterocycle_formation_reactions, snar_diversification_reactions, reaction_sequence, findings_json

        if node["type"] == "reaction" and "metadata" in node and "mapped_reaction_smiles" in node["metadata"]:
            rsmi = node["metadata"]["mapped_reaction_smiles"]
            reactants = rsmi.split(">")[0].split(".")
            product = rsmi.split(">")[-1]

            print(f"Examining reaction at depth {depth}: {rsmi}")

            # Check for acid chloride formation
            if any(checker.check_fg("Carboxylic acid", r) for r in reactants):
                if "Carboxylic acid" not in findings_json["atomic_checks"]["functional_groups"]:
                    findings_json["atomic_checks"]["functional_groups"].append("Carboxylic acid")
            if checker.check_fg("Acyl halide", product):
                if "Acyl halide" not in findings_json["atomic_checks"]["functional_groups"]:
                    findings_json["atomic_checks"]["functional_groups"].append("Acyl halide")

            if any(checker.check_fg("Carboxylic acid", r) for r in reactants) and checker.check_fg(
                "Acyl halide", product
            ):
                print("✓ Detected acid chloride formation")
                acid_chloride_reactions.append((depth, rsmi))

            # Check for amide formation
            has_acid_chloride = any(checker.check_fg("Acyl halide", r) for r in reactants)
            if has_acid_chloride and "Acyl halide" not in findings_json["atomic_checks"]["functional_groups"]:
                findings_json["atomic_checks"]["functional_groups"].append("Acyl halide")

            has_amine = False
            for r in reactants:
                if checker.check_fg("Aniline", r):
                    has_amine = True
                    if "Aniline" not in findings_json["atomic_checks"]["functional_groups"]:
                        findings_json["atomic_checks"]["functional_groups"].append("Aniline")
                if checker.check_fg("Primary amine", r):
                    has_amine = True
                    if "Primary amine" not in findings_json["atomic_checks"]["functional_groups"]:
                        findings_json["atomic_checks"]["functional_groups"].append("Primary amine")
                if checker.check_fg("Secondary amine", r):
                    has_amine = True
                    if "Secondary amine" not in findings_json["atomic_checks"]["functional_groups"]:
                        findings_json["atomic_checks"]["functional_groups"].append("Secondary amine")

            has_amide = False
            if checker.check_fg("Primary amide", product):
                has_amide = True
                if "Primary amide" not in findings_json["atomic_checks"]["functional_groups"]:
                    findings_json["atomic_checks"]["functional_groups"].append("Primary amide")
            if checker.check_fg("Secondary amide", product):
                has_amide = True
                if "Secondary amide" not in findings_json["atomic_checks"]["functional_groups"]:
                    findings_json["atomic_checks"]["functional_groups"].append("Secondary amide")
            if checker.check_fg("Tertiary amide", product):
                has_amide = True
                if "Tertiary amide" not in findings_json["atomic_checks"]["functional_groups"]:
                    findings_json["atomic_checks"]["functional_groups"].append("Tertiary amide")

            if has_acid_chloride and has_amine and has_amide:
                print("✓ Detected amide formation")
                amide_formation_reactions.append((depth, rsmi))
            
            # Specific reaction checks for amide formation
            amide_rxn_names = [
                "Acyl chloride with primary amine to amide (Schotten-Baumann)",
                "Acyl chloride with secondary amine to amide",
                "Schotten-Baumann_amide",
                "Acylation of Nitrogen Nucleophiles by Acyl/Thioacyl/Carbamoyl Halides and Analogs_N",
            ]
            for r_name in amide_rxn_names:
                if checker.check_reaction(r_name, rsmi):
                    print(f"✓ Detected amide formation reaction: {r_name}")
                    amide_formation_reactions.append((depth, rsmi))
                    if r_name not in findings_json["atomic_checks"]["named_reactions"]:
                        findings_json["atomic_checks"]["named_reactions"].append(r_name)
                    break # Only need to find one

            # Check for heterocycle formation
            product_has_heterocycle = False
            reactants_have_heterocycle = False

            for ring in HETEROCYCLES_OF_INTEREST:
                if checker.check_ring(ring, product):
                    product_has_heterocycle = True
                    if ring not in findings_json["atomic_checks"]["ring_systems"]:
                        findings_json["atomic_checks"]["ring_systems"].append(ring)
                if any(checker.check_ring(ring, r) for r in reactants):
                    reactants_have_heterocycle = True
                    if ring not in findings_json["atomic_checks"]["ring_systems"]:
                        findings_json["atomic_checks"]["ring_systems"].append(ring)

            if product_has_heterocycle and not reactants_have_heterocycle:
                print("✓ Detected heterocycle formation")
                heterocycle_formation_reactions.append((depth, rsmi))
            
            # Specific reaction checks for heterocycle formation
            heterocycle_rxn_names = [
                "Benzoxazole formation from aldehyde",
                "Benzoxazole formation from acyl halide",
                "Benzoxazole formation from ester/carboxylic acid",
                "Benzimidazole formation from aldehyde",
                "Benzimidazole formation from acyl halide",
                "Benzimidazole formation from ester/carboxylic acid",
                "Benzothiazole formation from aldehyde",
                "Benzothiazole formation from acyl halide",
                "Benzothiazole formation from ester/carboxylic acid",
            ]
            for r_name in heterocycle_rxn_names:
                if checker.check_reaction(r_name, rsmi):
                    print(f"✓ Detected heterocycle formation reaction: {r_name}")
                    heterocycle_formation_reactions.append((depth, rsmi))
                    if r_name not in findings_json["atomic_checks"]["named_reactions"]:
                        findings_json["atomic_checks"]["named_reactions"].append(r_name)
                    break # Only need to find one

            # Check for SNAr diversification
            # First check if the product contains a heterocycle
            product_has_heterocycle_for_snar = False
            for ring in HETEROCYCLES_OF_INTEREST:
                if checker.check_ring(ring, product):
                    product_has_heterocycle_for_snar = True
                    if ring not in findings_json["atomic_checks"]["ring_systems"]:
                        findings_json["atomic_checks"]["ring_systems"].append(ring)
                    break

            if product_has_heterocycle_for_snar:
                # Then check for SNAr reactions
                snar_rxn_names = [
                    "N-arylation (Buchwald-Hartwig/Ullmann-Goldberg)",
                    "Buchwald-Hartwig/Ullmann-Goldberg/N-arylation primary amine",
                    "Buchwald-Hartwig/Ullmann-Goldberg/N-arylation secondary amine",
                    "Buchwald-Hartwig",
                    "Ullmann-Goldberg Substitution amine",
                    "Goldberg coupling",
                ]
                for r_name in snar_rxn_names:
                    if checker.check_reaction(r_name, rsmi):
                        print(f"✓ Detected SNAr diversification: {r_name}")
                        snar_diversification_reactions.append((depth, rsmi))
                        if r_name not in findings_json["atomic_checks"]["named_reactions"]:
                            findings_json["atomic_checks"]["named_reactions"].append(r_name)
                        break # Only need to find one

            # Add to reaction sequence
            reaction_type = None
            if acid_chloride_reactions and (depth, rsmi) == acid_chloride_reactions[-1]:
                reaction_type = "acid_chloride"
            elif amide_formation_reactions and (depth, rsmi) == amide_formation_reactions[-1]:
                reaction_type = "amide_formation"
            elif (
                heterocycle_formation_reactions
                and (depth, rsmi) == heterocycle_formation_reactions[-1]
            ):
                reaction_type = "heterocycle_formation"
            elif (
                snar_diversification_reactions
                and (depth, rsmi) == snar_diversification_reactions[-1]
            ):
                reaction_type = "snar_diversification"

            if reaction_type:
                reaction_sequence.append((depth, reaction_type))

        # Continue traversing
        for child in node.get("children", []):
            # New logic for depth calculation
            new_depth = depth
            if node["type"] != "reaction": # If current node is chemical, depth increases
                new_depth = depth + 1
            # If current node is reaction, depth remains the same
            dfs_traverse(child, new_depth)

    # Start traversal
    dfs_traverse(route)

    # Sort reactions by depth to determine sequence
    reaction_sequence.sort(key=lambda x: x[0])

    # Print summary of detected steps
    print("\nSynthesis Strategy Detection Summary:")
    print(f"Acid chloride formation: {'✓' if acid_chloride_reactions else '✗'}")
    print(f"Amide formation: {'✓' if amide_formation_reactions else '✗'}")
    print(f"Heterocycle formation: {'✓' if heterocycle_formation_reactions else '✗'}")
    print(f"SNAr diversification: {'✓' if snar_diversification_reactions else '✗'}")

    if reaction_sequence:
        print("\nReaction Sequence (by depth):")
        for depth, reaction_type in reaction_sequence:
            print(f"Depth {depth}: {reaction_type}")

    # Check if all steps are present
    has_all_steps = bool(
        acid_chloride_reactions
        and amide_formation_reactions
        and heterocycle_formation_reactions
        and snar_diversification_reactions
    )

    # Check if the sequence is correct (in retrosynthetic order)
    correct_sequence = False
    if has_all_steps and len(reaction_sequence) >= 4:
        # Extract just the reaction types in sequence
        types_in_sequence = [rt for _, rt in reaction_sequence]

        # In retrosynthetic traversal, SNAr should have lowest depth (earliest in traversal)
        # and acid chloride should have highest depth (latest in traversal)
        snar_idx = (
            types_in_sequence.index("snar_diversification")
            if "snar_diversification" in types_in_sequence
            else -1
        )
        heterocycle_idx = (
            types_in_sequence.index("heterocycle_formation")
            if "heterocycle_formation" in types_in_sequence
            else -1
        )
        amide_idx = (
            types_in_sequence.index("amide_formation")
            if "amide_formation" in types_in_sequence
            else -1
        )
        acid_chloride_idx = (
            types_in_sequence.index("acid_chloride") if "acid_chloride" in types_in_sequence else -1
        )

        # Check if indices are in the correct order
        if snar_idx >= 0 and heterocycle_idx >= 0 and amide_idx >= 0 and acid_chloride_idx >= 0:
            # In retrosynthetic traversal, SNAr should have lowest depth (earliest in traversal)
            # and acid chloride should have highest depth (latest in traversal)
            correct_sequence = snar_idx <= heterocycle_idx <= amide_idx <= acid_chloride_idx

    print(f"\nAll steps present: {'✓' if has_all_steps else '✗'}")
    print(f"Correct sequence: {'✓' if correct_sequence else '✗'}")

    # Populate structural constraints in findings_json
    if has_all_steps:
        findings_json["structural_constraints"].append({
            "type": "co-occurrence",
            "details": {
                "targets": [
                    "acid_chloride_formation",
                    "amide_formation",
                    "heterocycle_formation",
                    "snar_diversification"
                ],
                "min_occurrences": 1,
                "description": "The synthesis must contain all four key reaction types: acid chloride formation, amide formation, heterocycle formation, and SNAr diversification."
            }
        })
    
    if correct_sequence:
        findings_json["structural_constraints"].append({
            "type": "sequence",
            "details": {
                "ordered_steps": [
                    "snar_diversification",
                    "heterocycle_formation",
                    "amide_formation",
                    "acid_chloride_formation"
                ],
                "description": "The key steps must occur in a specific order during the retrosynthetic traversal (ordered by increasing depth). This corresponds to a forward synthesis sequence of Acid Chloride Formation -> Amide Formation -> Heterocycle Formation -> SNAr Diversification."
            }
        })

    # Return True if all key steps are detected in the correct sequence
    result = has_all_steps and correct_sequence
    return result, findings_json
