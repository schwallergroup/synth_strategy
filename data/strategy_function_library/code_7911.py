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


N_HETEROCYCLES_OF_INTEREST = [
    "pyridine", "pyrazole", "imidazole", "oxazole", "thiazole",
    "pyrimidine", "pyrazine", "pyridazine", "triazole", "tetrazole",
    "indole", "quinoline", "isoquinoline", "purine", "benzimidazole",
    "benzoxazole", "benzothiazole", "indazole", "benzotriazole"
]

AMIDE_FORMATION_REACTIONS = [
    "Ester with primary amine to amide",
    "Ester with secondary amine to amide",
    "Aminolysis of esters",
    "Carboxylic acid with primary amine to amide",
    "Acyl chloride with primary amine to amide (Schotten-Baumann)",
    "Acyl chloride with secondary amine to amide",
    "Schotten-Baumann_amide",
    "Acylation of primary amines",
    "Acylation of secondary amines",
    "Acylation of Nitrogen Nucleophiles by Carboxylic Acids",
]

SNAR_N_ARYLATION_REACTIONS = [
    "heteroaromatic_nuc_sub",
    "nucl_sub_aromatic_ortho_nitro",
    "nucl_sub_aromatic_para_nitro",
    "N-arylation (Buchwald-Hartwig/Ullmann-Goldberg)",
    "Buchwald-Hartwig/Ullmann-Goldberg/N-arylation primary amine",
    "Buchwald-Hartwig/Ullmann-Goldberg/N-arylation secondary amine",
    "Buchwald-Hartwig",
    "N-arylation_heterocycles",
    "Goldberg coupling",
    "Ullmann-Goldberg Substitution amine",
]

def main(route) -> Tuple[bool, Dict]:
    """
    This function detects a strategy involving early amide formation from ester and amine,
    followed by late-stage SNAr reaction to incorporate a nitrogen heterocycle.
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
    has_early_amide_formation = False
    has_late_stage_snar = False

    # Store all reactions in order of traversal
    all_reactions = []

    # Calculate depth of a node in the synthesis tree
    def calculate_depth(node, root):
        if node == root:
            return 0

        def find_depth_recursive(current, target, current_depth=0):
            if current == target:
                return current_depth

            for child in current.get("children", []):
                result = find_depth_recursive(child, target, current_depth + 1)
                if result is not None:
                    return result
            return None

        return find_depth_recursive(root, node)

    def dfs_traverse(node, depth=0):
        nonlocal has_early_amide_formation, has_late_stage_snar, findings_json

        if node["type"] == "reaction":
            # For this implementation, we'll use the passed depth
            current_depth = depth

            rsmi = node.get("metadata", {}).get("rsmi", "")

            if rsmi:
                all_reactions.append((current_depth, rsmi))
                reactants = rsmi.split(">")[0].split(".")
                product = rsmi.split(">")[-1]

                print(f"Analyzing reaction at depth {current_depth}: {rsmi}")

                # Check for amide formation
                is_amide_formation = False
                for rxn in AMIDE_FORMATION_REACTIONS:
                    if checker.check_reaction(rxn, rsmi):
                        is_amide_formation = True
                        if rxn not in findings_json["atomic_checks"]["named_reactions"]:
                            findings_json["atomic_checks"]["named_reactions"].append(rxn)

                # Check if product has amide
                has_amide_product = False
                for fg in ["Secondary amide", "Primary amide", "Tertiary amide"]:
                    if checker.check_fg(fg, product):
                        has_amide_product = True
                        if fg not in findings_json["atomic_checks"]["functional_groups"]:
                            findings_json["atomic_checks"]["functional_groups"].append(fg)

                # Check if reactants include ester/acid and amine
                has_ester_or_acid = False
                for r in reactants:
                    for fg in ["Ester", "Carboxylic acid", "Acyl halide"]:
                        if checker.check_fg(fg, r):
                            has_ester_or_acid = True
                            if fg not in findings_json["atomic_checks"]["functional_groups"]:
                                findings_json["atomic_checks"]["functional_groups"].append(fg)

                has_amine = False
                for r in reactants:
                    for fg in ["Primary amine", "Secondary amine", "Aniline"]:
                        if checker.check_fg(fg, r):
                            has_amine = True
                            if fg not in findings_json["atomic_checks"]["functional_groups"]:
                                findings_json["atomic_checks"]["functional_groups"].append(fg)

                # Check if reactants don't have amide but product does (amide formation)
                reactants_have_amide = False
                for r in reactants:
                    for fg in ["Primary amide", "Secondary amide", "Tertiary amide"]:
                        if checker.check_fg(fg, r):
                            reactants_have_amide = True
                            # No need to add to findings_json here, as it's a negative condition

                if is_amide_formation or (
                    has_amide_product
                    and has_ester_or_acid
                    and has_amine
                    and not reactants_have_amide
                ):
                    print(f"✓ Detected amide formation at depth {current_depth}")
                    if "amide_formation" not in findings_json["atomic_checks"]["named_reactions"]:
                        findings_json["atomic_checks"]["named_reactions"].append("amide_formation")
                    # Mark as early if depth >= 4
                    if current_depth >= 4:
                        print(f"✓ This is considered early-stage amide formation")
                        has_early_amide_formation = True
                        if {"type": "positional", "details": {"target": "amide_formation", "position_check": "depth >= 4", "description": "An amide formation reaction (either from a named list or detected by functional group changes) must occur in an early stage of the synthesis, defined as a depth of 4 or greater."}} not in findings_json["structural_constraints"]:
                            findings_json["structural_constraints"].append({"type": "positional", "details": {"target": "amide_formation", "position_check": "depth >= 4", "description": "An amide formation reaction (either from a named list or detected by functional group changes) must occur in an early stage of the synthesis, defined as a depth of 4 or greater."}})

                # Check for SNAr or N-arylation reactions
                is_snar_or_arylation = False
                for rxn in SNAR_N_ARYLATION_REACTIONS:
                    if checker.check_reaction(rxn, rsmi):
                        is_snar_or_arylation = True
                        if rxn not in findings_json["atomic_checks"]["named_reactions"]:
                            findings_json["atomic_checks"]["named_reactions"].append(rxn)
                
                if is_snar_or_arylation and "SNAr_N_arylation_reaction" not in findings_json["atomic_checks"]["named_reactions"]:
                    findings_json["atomic_checks"]["named_reactions"].append("SNAr_N_arylation_reaction")

                # Check for nitrogen heterocycles in product
                product_n_heterocycles = []
                for ring in N_HETEROCYCLES_OF_INTEREST:
                    if checker.check_ring(ring, product):
                        product_n_heterocycles.append(ring)
                        if ring not in findings_json["atomic_checks"]["ring_systems"]:
                            findings_json["atomic_checks"]["ring_systems"].append(ring)

                if is_snar_or_arylation and product_n_heterocycles:
                    print(
                        f"✓ Detected SNAr or N-arylation at depth {current_depth} with nitrogen heterocycle: {product_n_heterocycles}"
                    )
                    if {"type": "co-occurrence", "details": {"targets": ["SNAr_N_arylation_reaction", "N_heterocycle_in_product"], "scope": "single_reaction", "description": "A target late-stage reaction is identified if it is a named SNAr/N-arylation reaction AND its product contains a specified nitrogen heterocycle."}} not in findings_json["structural_constraints"]:
                        findings_json["structural_constraints"].append({"type": "co-occurrence", "details": {"targets": ["SNAr_N_arylation_reaction", "N_heterocycle_in_product"], "scope": "single_reaction", "description": "A target late-stage reaction is identified if it is a named SNAr/N-arylation reaction AND its product contains a specified nitrogen heterocycle."}})
                    if "N_heterocycle_in_product" not in findings_json["atomic_checks"]["named_reactions"]:
                        findings_json["atomic_checks"]["named_reactions"].append("N_heterocycle_in_product")

                    # Mark as late-stage if depth <= 2
                    if current_depth <= 2:
                        print(f"✓ This is considered late-stage SNAr")
                        has_late_stage_snar = True
                        if {"type": "positional", "details": {"target": "SNAr_N_arylation_with_N_heterocycle", "position_check": "depth <= 2", "description": "An SNAr/N-arylation reaction that incorporates a nitrogen heterocycle must occur in a late stage of the synthesis, defined as a depth of 2 or less."}} not in findings_json["structural_constraints"]:
                            findings_json["structural_constraints"].append({"type": "positional", "details": {"target": "SNAr_N_arylation_with_N_heterocycle", "position_check": "depth <= 2", "description": "An SNAr/N-arylation reaction that incorporates a nitrogen heterocycle must occur in a late stage of the synthesis, defined as a depth of 2 or less."}})

        # Determine the depth for the recursive call
        next_depth = depth
        if node["type"] == "chemical": # Depth increases when going from chemical to reaction
            next_depth = depth + 1

        # Traverse children
        for child in node.get("children", []):
            dfs_traverse(child, next_depth)

    # Start traversal
    dfs_traverse(route, depth=0)

    # If we have reactions but couldn't determine early/late stages properly,
    # try to infer the strategy from the sequence of reactions
    if all_reactions and not (has_early_amide_formation and has_late_stage_snar):
        # Sort reactions by depth/position
        all_reactions.sort(key=lambda x: x[0])

        amide_indices = []
        snar_indices = []

        for i, (depth, rsmi) in enumerate(all_reactions):
            reactants = rsmi.split(">")[0].split(".")
            product = rsmi.split(">")[-1]

            # Check for amide formation
            is_amide_formation = any(
                checker.check_reaction(rxn, rsmi) for rxn in AMIDE_FORMATION_REACTIONS
            )

            # Check if product has amide
            has_amide_product = (
                checker.check_fg("Secondary amide", product)
                or checker.check_fg("Primary amide", product)
                or checker.check_fg("Tertiary amide", product)
            )

            # Check if reactants include ester/acid and amine
            has_ester_or_acid = any(
                checker.check_fg(fg, r)
                for r in reactants
                for fg in ["Ester", "Carboxylic acid", "Acyl halide"]
            )
            has_amine = any(
                checker.check_fg(fg, r)
                for r in reactants
                for fg in ["Primary amine", "Secondary amine", "Aniline"]
            )

            # Check if reactants don't have amide but product does (amide formation)
            reactants_have_amide = any(
                checker.check_fg(fg, r)
                for r in reactants
                for fg in ["Primary amide", "Secondary amide", "Tertiary amide"]
            )

            if is_amide_formation or (
                has_amide_product and has_ester_or_acid and has_amine and not reactants_have_amide
            ):
                amide_indices.append(i)

            # Check for SNAr or N-arylation
            is_snar_or_arylation = any(checker.check_reaction(rxn, rsmi) for rxn in SNAR_N_ARYLATION_REACTIONS)

            # Check for nitrogen heterocycles in product
            product_n_heterocycles = [
                ring for ring in N_HETEROCYCLES_OF_INTEREST if checker.check_ring(ring, product)
            ]

            if is_snar_or_arylation and product_n_heterocycles:
                snar_indices.append(i)

        # Check if we have both reaction types and amide formation comes before SNAr
        if amide_indices and snar_indices:
            min_amide_idx = min(amide_indices)
            max_snar_idx = max(snar_indices)

            if min_amide_idx < max_snar_idx:
                print(
                    f"Strategy detected from reaction sequence: amide formation at position {min_amide_idx} followed by SNAr at position {max_snar_idx}"
                )
                has_early_amide_formation = True
                has_late_stage_snar = True
                if {"type": "sequence", "details": {"ordered_events": ["amide_formation", "SNAr_N_arylation_with_N_heterocycle"], "description": "As a fallback condition, the strategy is detected if any amide formation occurs sequentially before any SNAr/N-arylation reaction that incorporates a nitrogen heterocycle, irrespective of depth."}} not in findings_json["structural_constraints"]:
                    findings_json["structural_constraints"].append({"type": "sequence", "details": {"ordered_events": ["amide_formation", "SNAr_N_arylation_with_N_heterocycle"], "description": "As a fallback condition, the strategy is detected if any amide formation occurs sequentially before any SNAr/N-arylation reaction that incorporates a nitrogen heterocycle, irrespective of depth."}})

    # Return True if both key elements of the strategy are present
    result = has_early_amide_formation and has_late_stage_snar
    print(
        f"Strategy detected: early amide formation: {has_early_amide_formation}, late-stage SNAr: {has_late_stage_snar}"
    )
    return result, findings_json
