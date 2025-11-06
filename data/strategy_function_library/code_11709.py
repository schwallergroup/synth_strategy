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


# Refactored lists for enumeration
SUZUKI_REACTIONS = [
    "Suzuki coupling with boronic acids",
    "Suzuki coupling with boronic acids OTf",
    "Suzuki coupling with sulfonic esters",
    "Suzuki coupling with boronic esters",
    "Suzuki coupling with boronic esters OTf",
    "{Suzuki}",
]
OTHER_CROSS_COUPLINGS = [
    "Negishi coupling",
    "Stille reaction_aryl",
    "Hiyama-Denmark Coupling",
    "Kumada cross-coupling",
    "Ullmann condensation",
    "Buchwald-Hartwig",
    "{Buchwald-Hartwig}",
]
SULFONYL_FGS = ["Sulfonic acid", "Sulfonate", "Triflate", "Mesylate", "Tosylate"]
SULFONATION_REACTIONS = [
    "Formation of Sulfonic Esters",
    "Formation of Sulfonic Esters on TMS protected alcohol",
    "Alcohol to triflate conversion",
]
HALIDE_FGS = ["Aromatic halide", "Primary halide", "Secondary halide", "Tertiary halide"]
HALOGENATION_REACTIONS = [
    "Aromatic fluorination",
    "Aromatic chlorination",
    "Aromatic bromination",
    "Aromatic iodination",
    "Chlorination",
    "Fluorination",
    "Iodination",
    "Bromination",
    "Wohl-Ziegler bromination benzyl primary",
    "Wohl-Ziegler bromination benzyl secondary",
    "Wohl-Ziegler bromination benzyl tertiary",
    "Wohl-Ziegler bromination allyl primary",
    "Alcohol to chloride_sulfonyl chloride",
    "Alcohol to chloride_SOCl2",
    "Primary amine to chloride",
]

def main(route) -> Tuple[bool, Dict]:
    """
    Detects a strategy of early-stage handle installation (sulfonation or halogenation)
    followed by late-stage cross-coupling. Cross-coupling is identified by reaction name
    (see SUZUKI_REACTIONS, OTHER_CROSS_COUPLINGS) or by the presence of specific reactants
    (e.g., boronic acids and aryl halides). Handle installation is identified by reaction
    name (see SULFONATION_REACTIONS, HALOGENATION_REACTIONS) or by the net formation
    of a sulfonyl or halide functional group.
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

    # Track if we found the key features
    found_cross_coupling = False
    found_sulfonation = False
    found_halogenation = False

    # Track the depth at which each transformation occurs
    cross_coupling_depth = None
    sulfonation_depth = None
    halogenation_depths = []

    def dfs_traverse(node, depth=0):
        nonlocal found_cross_coupling, found_sulfonation, found_halogenation
        nonlocal cross_coupling_depth, sulfonation_depth, halogenation_depths, findings_json

        if node["type"] == "reaction":
            try:
                # Extract reaction SMILES
                rsmi = node["metadata"]["mapped_reaction_smiles"]
                reactants_part = rsmi.split(">")[0]
                products_part = rsmi.split(">")[-1]
                reactants = reactants_part.split(".")
                product = products_part  # Usually a single product

                # Check for Suzuki coupling by reaction type or by reactant functional groups
                for rxn_type in SUZUKI_REACTIONS:
                    if checker.check_reaction(rxn_type, rsmi):
                        found_cross_coupling = True
                        cross_coupling_depth = depth
                        findings_json["atomic_checks"]["named_reactions"].append(rxn_type)
                        print(f"Found Suzuki cross-coupling at depth {depth}")
                        break
                
                boronic_acid_found_in_reactants = False
                for r in reactants:
                    if checker.check_fg("Boronic acid", r):
                        findings_json["atomic_checks"]["functional_groups"].append("Boronic acid")
                        boronic_acid_found_in_reactants = True
                    if checker.check_fg("Boronic ester", r):
                        findings_json["atomic_checks"]["functional_groups"].append("Boronic ester")
                        boronic_acid_found_in_reactants = True

                aryl_halide_or_triflate_found_in_reactants = False
                for r in reactants:
                    if checker.check_fg("Aromatic halide", r):
                        findings_json["atomic_checks"]["functional_groups"].append("Aromatic halide")
                        aryl_halide_or_triflate_found_in_reactants = True
                    if checker.check_fg("Triflate", r):
                        findings_json["atomic_checks"]["functional_groups"].append("Triflate")
                        aryl_halide_or_triflate_found_in_reactants = True

                if not found_cross_coupling and boronic_acid_found_in_reactants and aryl_halide_or_triflate_found_in_reactants:
                    found_cross_coupling = True
                    cross_coupling_depth = depth
                    print(f"Found Suzuki cross-coupling (by reactants) at depth {depth}")

                # Check for other cross-coupling reactions that might form biaryls
                if not found_cross_coupling:
                    for rxn_type in OTHER_CROSS_COUPLINGS:
                        if checker.check_reaction(rxn_type, rsmi):
                            found_cross_coupling = True
                            cross_coupling_depth = depth
                            findings_json["atomic_checks"]["named_reactions"].append(rxn_type)
                            print(f"Found alternative cross-coupling at depth {depth}")
                            break

                # Check for sulfonation (introduction of sulfonic acid group or sulfonate derivatives)
                product_has_sulfonic = False
                for fg in SULFONYL_FGS:
                    if checker.check_fg(fg, product):
                        product_has_sulfonic = True
                        if fg not in findings_json["atomic_checks"]["functional_groups"]:
                            findings_json["atomic_checks"]["functional_groups"].append(fg)

                reactants_have_sulfonic = False
                for r in reactants:
                    for fg in SULFONYL_FGS:
                        if checker.check_fg(fg, r):
                            reactants_have_sulfonic = True
                            if fg not in findings_json["atomic_checks"]["functional_groups"]:
                                findings_json["atomic_checks"]["functional_groups"].append(fg)

                sulfonation_by_reaction = False
                for rxn_type in SULFONATION_REACTIONS:
                    if checker.check_reaction(rxn_type, rsmi):
                        sulfonation_by_reaction = True
                        findings_json["atomic_checks"]["named_reactions"].append(rxn_type)
                        break

                if (product_has_sulfonic and not reactants_have_sulfonic):
                    found_sulfonation = True
                    sulfonation_depth = depth
                    findings_json["atomic_checks"]["named_reactions"].append("sulfonation_by_net_change")
                    print(f"Found sulfonation at depth {depth} (by net change)")
                
                if sulfonation_by_reaction:
                    found_sulfonation = True
                    sulfonation_depth = depth
                    print(f"Found sulfonation at depth {depth} (by reaction)")

                # Check for halogenation (introduction of halides)
                product_has_halide = False
                for fg in HALIDE_FGS:
                    if checker.check_fg(fg, product):
                        product_has_halide = True
                        if fg not in findings_json["atomic_checks"]["functional_groups"]:
                            findings_json["atomic_checks"]["functional_groups"].append(fg)

                reactants_have_halide = False
                for r in reactants:
                    for fg in HALIDE_FGS:
                        if checker.check_fg(fg, r):
                            reactants_have_halide = True
                            if fg not in findings_json["atomic_checks"]["functional_groups"]:
                                findings_json["atomic_checks"]["functional_groups"].append(fg)

                halogenation_by_reaction = False
                for rxn_type in HALOGENATION_REACTIONS:
                    if checker.check_reaction(rxn_type, rsmi):
                        halogenation_by_reaction = True
                        findings_json["atomic_checks"]["named_reactions"].append(rxn_type)
                        break

                if (product_has_halide and not reactants_have_halide):
                    found_halogenation = True
                    halogenation_depths.append(depth)
                    findings_json["atomic_checks"]["named_reactions"].append("halogenation_by_net_change")
                    print(f"Found halogenation at depth {depth} (by net change)")
                
                if halogenation_by_reaction:
                    found_halogenation = True
                    halogenation_depths.append(depth)
                    print(f"Found halogenation at depth {depth} (by reaction)")

            except Exception as e:
                print(f"Error processing reaction node: {e}")

        # Traverse children
        for child in node.get("children", []):
            # New logic for depth calculation
            if node["type"] == "reaction":
                # Depth remains the same when traversing from a reaction node to a chemical node
                dfs_traverse(child, depth)
            else:
                # Depth increases when traversing from a chemical node to a reaction node
                dfs_traverse(child, depth + 1)

    # Start traversal
    dfs_traverse(route)

    # Check if the strategy is present
    strategy_present = found_cross_coupling and (found_sulfonation or found_halogenation)

    if strategy_present and cross_coupling_depth is not None:
        handle_installation_found = False
        if found_sulfonation and sulfonation_depth is not None:
            handle_installation_found = True
            if cross_coupling_depth < sulfonation_depth:
                strategy_present = False # Cross-coupling should be AFTER handle installation
            print(
                f"Strategy assessment: Cross-coupling at depth {cross_coupling_depth}, Sulfonation at depth {sulfonation_depth}"
            )

        if found_halogenation and halogenation_depths:
            handle_installation_found = True
            earliest_halogenation = min(halogenation_depths)
            if cross_coupling_depth < earliest_halogenation:
                strategy_present = False # Cross-coupling should be AFTER handle installation
            print(
                f"Strategy assessment: Cross-coupling at depth {cross_coupling_depth}, Earliest halogenation at depth {earliest_halogenation}"
            )
        
        if found_cross_coupling and handle_installation_found:
            findings_json["structural_constraints"].append({
                "type": "co-occurrence",
                "details": {
                    "targets": [
                        "cross_coupling",
                        "handle_installation"
                    ]
                }
            })
            # Re-evaluate sequence constraint based on the correct logic (handle_installation BEFORE cross_coupling)
            if (found_sulfonation and sulfonation_depth is not None and cross_coupling_depth > sulfonation_depth) or \
               (found_halogenation and halogenation_depths and cross_coupling_depth > min(halogenation_depths)):
                findings_json["structural_constraints"].append({
                    "type": "sequence",
                    "details": {
                        "before": "handle_installation",
                        "after": "cross_coupling"
                    }
                })
            else:
                strategy_present = False # Sequence constraint not met

    print(f"Final strategy assessment: {strategy_present}")
    return strategy_present, findings_json
