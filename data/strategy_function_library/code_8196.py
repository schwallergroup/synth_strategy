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


AROMATIC_HALOGENATION_REACTIONS = [
    "Aromatic fluorination",
    "Aromatic chlorination",
    "Aromatic bromination",
    "Aromatic iodination",
]

CROSS_COUPLING_REACTIONS = [
    "Suzuki coupling with boronic acids",
    "Suzuki coupling with boronic acids OTf",
    "Suzuki coupling with sulfonic esters",
    "Suzuki coupling with boronic esters OTf",
    "Suzuki coupling with boronic esters",
    "Negishi coupling",
    "Stille reaction_vinyl",
    "Stille reaction_aryl",
    "Stille reaction_benzyl",
    "Stille reaction_allyl",
    "Stille reaction_vinyl OTf",
    "Stille reaction_aryl OTf",
    "Stille reaction_benzyl OTf",
    "Stille reaction_allyl OTf",
    "Stille reaction_other",
    "Stille reaction_other OTf",
    "Hiyama-Denmark Coupling",
    "Kumada cross-coupling",
    "Aryllithium cross-coupling",
    "{Suzuki}",
    "Buchwald-Hartwig/Ullmann-Goldberg/N-arylation primary amine",
    "Buchwald-Hartwig/Ullmann-Goldberg/N-arylation secondary amine",
    "{Buchwald-Hartwig}",
    "N-arylation (Buchwald-Hartwig/Ullmann-Goldberg)",
    "Ullmann-Goldberg Substitution amine",
    "Goldberg coupling aryl amine-aryl chloride",
    "Goldberg coupling aryl amide-aryl chloride",
    "Goldberg coupling",
    "Ullmann condensation",
]

def main(route) -> Tuple[bool, Dict]:
    """
    Detects a late-stage aromatic cross-coupling strategy. This is defined as an early-stage (depth > 1) aromatic halogenation followed by a late-stage (depth <= 1) cross-coupling. Aromatic halogenation is identified using named reactions from the `AROMATIC_HALOGENATION_REACTIONS` list or by the net formation of an aromatic halide. Cross-coupling is identified using named reactions from the `CROSS_COUPLING_REACTIONS` list or by the consumption of an aromatic halide with a common organometallic coupling partner.
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

    # Track if we found a cross-coupling reaction at depth 0 or 1 (late stage)
    found_cross_coupling = False
    # Track if we found a halogenation reaction in an earlier step
    found_halogenation = False
    # Track the depth of halogenation
    halogenation_depth = -1
    # Track the depth of cross-coupling
    cross_coupling_depth = -1

    def dfs_traverse(node, depth=0):
        nonlocal found_cross_coupling, found_halogenation, halogenation_depth, cross_coupling_depth, findings_json

        if node["type"] == "reaction":
            # Extract reactants and product
            rsmi = node["metadata"].get("rsmi", "")
            if not rsmi:
                return

            reactants = rsmi.split(">")[0].split(".")
            product = rsmi.split(">")[-1]

            print(f"Analyzing reaction at depth {depth}: {rsmi}")

            # Check if this is a cross-coupling reaction (depth 0 or 1 = late stage)
            if depth <= 1:
                # First check for named cross-coupling reactions
                for reaction_type in CROSS_COUPLING_REACTIONS:
                    if checker.check_reaction(reaction_type, rsmi):
                        print(f"Detected {reaction_type} reaction at depth {depth}")
                        found_cross_coupling = True
                        cross_coupling_depth = depth
                        if reaction_type not in findings_json["atomic_checks"]["named_reactions"]:
                            findings_json["atomic_checks"]["named_reactions"].append(reaction_type)
                        break

                # If no named reaction was found, check for general cross-coupling pattern
                if not found_cross_coupling:
                    # Check if any reactant has an aromatic halide that's not in the product
                    has_aromatic_halide_in_reactants = False
                    for reactant in reactants:
                        if checker.check_fg("Aromatic halide", reactant):
                            has_aromatic_halide_in_reactants = True
                            if "Aromatic halide" not in findings_json["atomic_checks"]["functional_groups"]:
                                findings_json["atomic_checks"]["functional_groups"].append("Aromatic halide")
                            print(f"Found aromatic halide in reactant: {reactant}")
                            break

                    has_aromatic_halide_in_product = checker.check_fg("Aromatic halide", product)
                    if has_aromatic_halide_in_product and "Aromatic halide" not in findings_json["atomic_checks"]["functional_groups"]:
                        findings_json["atomic_checks"]["functional_groups"].append("Aromatic halide")

                    # If aromatic halide is consumed in the reaction, it might be a cross-coupling
                    if has_aromatic_halide_in_reactants and not has_aromatic_halide_in_product:
                        # Check for common coupling partners
                        has_coupling_partner = False
                        for reactant in reactants:
                            if (
                                checker.check_fg("Boronic acid", reactant)
                                or checker.check_fg("Boronic ester", reactant)
                                or checker.check_fg("Magnesium halide", reactant)
                                or checker.check_fg("Zinc halide", reactant)
                                or "[Sn]" in reactant
                                or checker.check_fg("Tin", reactant)
                            ):
                                has_coupling_partner = True
                                if checker.check_fg("Boronic acid", reactant) and "Boronic acid" not in findings_json["atomic_checks"]["functional_groups"]:
                                    findings_json["atomic_checks"]["functional_groups"].append("Boronic acid")
                                if checker.check_fg("Boronic ester", reactant) and "Boronic ester" not in findings_json["atomic_checks"]["functional_groups"]:
                                    findings_json["atomic_checks"]["functional_groups"].append("Boronic ester")
                                if checker.check_fg("Magnesium halide", reactant) and "Magnesium halide" not in findings_json["atomic_checks"]["functional_groups"]:
                                    findings_json["atomic_checks"]["functional_groups"].append("Magnesium halide")
                                if checker.check_fg("Zinc halide", reactant) and "Zinc halide" not in findings_json["atomic_checks"]["functional_groups"]:
                                    findings_json["atomic_checks"]["functional_groups"].append("Zinc halide")
                                if checker.check_fg("Tin", reactant) and "Tin" not in findings_json["atomic_checks"]["functional_groups"]:
                                    findings_json["atomic_checks"]["functional_groups"].append("Tin")
                                print(f"Found coupling partner: {reactant}")
                                break

                        if has_coupling_partner:
                            found_cross_coupling = True
                            cross_coupling_depth = depth
                            if "generic_cross_coupling" not in findings_json["atomic_checks"]["named_reactions"]:
                                findings_json["atomic_checks"]["named_reactions"].append("generic_cross_coupling")
                            print(f"Found general cross-coupling pattern at depth {depth}")

            # Check if this is a halogenation reaction
            elif depth > 1:
                # Check for halogenation reactions
                for reaction_type in AROMATIC_HALOGENATION_REACTIONS:
                    if checker.check_reaction(reaction_type, rsmi):
                        found_halogenation = True
                        halogenation_depth = depth
                        if reaction_type not in findings_json["atomic_checks"]["named_reactions"]:
                            findings_json["atomic_checks"]["named_reactions"].append(reaction_type)
                        print(f"Found {reaction_type} reaction at depth {depth}")
                        break

                # If no specific halogenation reaction was found, check for appearance of aromatic halide
                if not found_halogenation:
                    has_halogen_in_product = checker.check_fg("Aromatic halide", product)
                    if has_halogen_in_product and "Aromatic halide" not in findings_json["atomic_checks"]["functional_groups"]:
                        findings_json["atomic_checks"]["functional_groups"].append("Aromatic halide")

                    has_halogen_in_reactants = False
                    for reactant in reactants:
                        if checker.check_fg("Aromatic halide", reactant):
                            has_halogen_in_reactants = True
                            if "Aromatic halide" not in findings_json["atomic_checks"]["functional_groups"]:
                                findings_json["atomic_checks"]["functional_groups"].append("Aromatic halide")
                            break

                    if has_halogen_in_product and not has_halogen_in_reactants:
                        found_halogenation = True
                        halogenation_depth = depth
                        if "aromatic_halogen_formation" not in findings_json["atomic_checks"]["named_reactions"]:
                            findings_json["atomic_checks"]["named_reactions"].append("aromatic_halogen_formation")
                        print(f"Found halogenation (functional group change) at depth {depth}")

        # Traverse children
        for child in node.get("children", []):
            # New depth calculation logic
            new_depth = depth
            if node["type"] != "reaction": # If current node is chemical, depth increases
                new_depth = depth + 1
            # If current node is reaction, depth remains the same
            dfs_traverse(child, new_depth)

    # Start traversal from the root
    dfs_traverse(route)

    # The strategy is present if we found both a cross-coupling at the end
    # and a halogenation step earlier in the synthesis
    strategy_present = found_cross_coupling and found_halogenation

    if strategy_present:
        print(
            f"Late-stage aromatic cross-coupling strategy detected: Halogenation at depth {halogenation_depth} followed by coupling at depth {cross_coupling_depth}"
        )
        # Add structural constraints if the strategy is detected
        findings_json["structural_constraints"].append({
            "type": "co-occurrence",
            "details": {
                "targets": [
                    "any_aromatic_halogenation",
                    "any_cross_coupling"
                ],
                "description": "The strategy requires both an aromatic halogenation event and a cross-coupling event to occur in the route."
            }
        })
        findings_json["structural_constraints"].append({
            "type": "positional",
            "details": {
                "target": "any_aromatic_halogenation",
                "position": "early_stage",
                "description": "An aromatic halogenation event (either a named reaction or net formation of an aromatic halide) must occur at a depth > 1."
            }
        })
        findings_json["structural_constraints"].append({
            "type": "positional",
            "details": {
                "target": "any_cross_coupling",
                "position": "late_stage",
                "description": "A cross-coupling event (either a named reaction or consumption of an aromatic halide with a coupling partner) must occur at a depth <= 1."
            }
        })
    else:
        print("Late-stage aromatic cross-coupling strategy not detected")
        if found_halogenation and not found_cross_coupling:
            print("Found halogenation but no cross-coupling at late stage")
            # If halogenation was found, but not cross-coupling, still record the halogenation positional constraint
            if halogenation_depth > 1:
                findings_json["structural_constraints"].append({
                    "type": "positional",
                    "details": {
                        "target": "any_aromatic_halogenation",
                        "position": "early_stage",
                        "description": "An aromatic halogenation event (either a named reaction or net formation of an aromatic halide) must occur at a depth > 1."
                    }
                })
        elif found_cross_coupling and not found_halogenation:
            print("Found cross-coupling at late stage but no halogenation in earlier steps")
            # If cross-coupling was found, but not halogenation, still record the cross-coupling positional constraint
            if cross_coupling_depth <= 1:
                findings_json["structural_constraints"].append({
                    "type": "positional",
                    "details": {
                        "target": "any_cross_coupling",
                        "position": "late_stage",
                        "description": "A cross-coupling event (either a named reaction or consumption of an aromatic halide with a coupling partner) must occur at a depth <= 1."
                    }
                })

    return strategy_present, findings_json
