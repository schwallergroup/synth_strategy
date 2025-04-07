#!/bin/python

"""LM-defined function for strategy description."""

import copy
import re
from collections import deque

import rdkit
import rdkit.Chem as Chem
from rdkit import Chem
from rdkit.Chem import (
    AllChem,
    Descriptors,
    Lipinski,
    rdChemReactions,
    rdFMCS,
    rdMolDescriptors,
    rdmolops,
)
from rdkit.Chem.Scaffolds import MurckoScaffold

from steerable_retro.utils import check, fuzzy_dict
from steerable_retro.utils.check import Check

fg_args = {
    "file_path": "/home/dparm/steerable_retro/data/patterns/functional_groups.json",
    "value_field": "pattern",
    "key_field": "name",
}
reaction_class_args = {
    "file_path": "/home/dparm/steerable_retro/data/patterns/smirks.json",
    "value_field": "smirks",
    "key_field": "name",
}
ring_smiles_args = {
    "file_path": "/home/dparm/steerable_retro/data/patterns/chemical_rings_smiles.json",
    "value_field": "smiles",
    "key_field": "name",
}
functional_groups = fuzzy_dict.FuzzyDict.from_json(**fg_args)
reaction_classes = fuzzy_dict.FuzzyDict.from_json(**reaction_class_args)
ring_smiles = fuzzy_dict.FuzzyDict.from_json(**ring_smiles_args)

checker = check.Check(
    fg_dict=functional_groups, reaction_dict=reaction_classes, ring_dict=ring_smiles
)


def main(route):
    """
    Detects if the synthesis route employs a late-stage aromatic cross-coupling strategy.
    This strategy involves halogenation of an aromatic ring in an early step,
    followed by a cross-coupling reaction in the final step.
    """
    # Track if we found a cross-coupling reaction at depth 0 or 1 (late stage)
    found_cross_coupling = False
    # Track if we found a halogenation reaction in an earlier step
    found_halogenation = False
    # Track the depth of halogenation
    halogenation_depth = -1
    # Track the depth of cross-coupling
    cross_coupling_depth = -1

    def dfs_traverse(node, depth=0):
        nonlocal found_cross_coupling, found_halogenation, halogenation_depth, cross_coupling_depth

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
                # Check for various cross-coupling reactions
                cross_coupling_reactions = [
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

                # First check for named cross-coupling reactions
                for reaction_type in cross_coupling_reactions:
                    if checker.check_reaction(reaction_type, rsmi):
                        print(f"Detected {reaction_type} reaction")

                        # For Stille reactions, check for organotin compounds
                        if "Stille" in reaction_type:
                            has_tin_compound = False
                            has_aromatic_halide = False

                            for reactant in reactants:
                                if "[Sn]" in reactant or checker.check_fg("Tin", reactant):
                                    has_tin_compound = True
                                    print(f"Found organotin compound: {reactant}")

                                if checker.check_fg("Aromatic halide", reactant):
                                    has_aromatic_halide = True
                                    print(f"Found aromatic halide: {reactant}")

                            if has_tin_compound and has_aromatic_halide:
                                found_cross_coupling = True
                                cross_coupling_depth = depth
                                print(f"Found Stille coupling at depth {depth}")
                                break
                        else:
                            # For other cross-coupling reactions, check for aromatic halide
                            for reactant in reactants:
                                if checker.check_fg("Aromatic halide", reactant):
                                    found_cross_coupling = True
                                    cross_coupling_depth = depth
                                    print(f"Found {reaction_type} reaction at depth {depth}")
                                    break

                        if found_cross_coupling:
                            break

                # If no named reaction was found, check for general cross-coupling pattern
                if not found_cross_coupling:
                    # Check if any reactant has an aromatic halide that's not in the product
                    has_aromatic_halide_in_reactants = False
                    for reactant in reactants:
                        if checker.check_fg("Aromatic halide", reactant):
                            has_aromatic_halide_in_reactants = True
                            print(f"Found aromatic halide in reactant: {reactant}")
                            break

                    has_aromatic_halide_in_product = checker.check_fg("Aromatic halide", product)

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
                                print(f"Found coupling partner: {reactant}")
                                break

                        if has_coupling_partner:
                            found_cross_coupling = True
                            cross_coupling_depth = depth
                            print(f"Found general cross-coupling pattern at depth {depth}")

            # Check if this is a halogenation reaction
            elif depth > 1:
                # Check for halogenation reactions
                halogenation_reactions = [
                    "Aromatic fluorination",
                    "Aromatic chlorination",
                    "Aromatic bromination",
                    "Aromatic iodination",
                ]

                for reaction_type in halogenation_reactions:
                    if checker.check_reaction(reaction_type, rsmi):
                        # Verify that product has aromatic halide but reactants don't
                        has_halogen_in_product = checker.check_fg("Aromatic halide", product)

                        has_halogen_in_reactants = False
                        for reactant in reactants:
                            if checker.check_fg("Aromatic halide", reactant):
                                has_halogen_in_reactants = True
                                break

                        if has_halogen_in_product and not has_halogen_in_reactants:
                            found_halogenation = True
                            halogenation_depth = depth
                            print(f"Found {reaction_type} reaction at depth {depth}")
                            break

                # If no specific halogenation reaction was found, check for appearance of aromatic halide
                if not found_halogenation:
                    has_halogen_in_product = checker.check_fg("Aromatic halide", product)

                    has_halogen_in_reactants = False
                    for reactant in reactants:
                        if checker.check_fg("Aromatic halide", reactant):
                            has_halogen_in_reactants = True
                            break

                    if has_halogen_in_product and not has_halogen_in_reactants:
                        found_halogenation = True
                        halogenation_depth = depth
                        print(f"Found halogenation (functional group change) at depth {depth}")

        # Traverse children
        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    # Start traversal from the root
    dfs_traverse(route)

    # The strategy is present if we found both a cross-coupling at the end
    # and a halogenation step earlier in the synthesis
    strategy_present = found_cross_coupling and found_halogenation

    if strategy_present:
        print(
            f"Late-stage aromatic cross-coupling strategy detected: Halogenation at depth {halogenation_depth} followed by coupling at depth {cross_coupling_depth}"
        )
    else:
        print("Late-stage aromatic cross-coupling strategy not detected")
        if found_halogenation and not found_cross_coupling:
            print("Found halogenation but no cross-coupling at late stage")
        elif found_cross_coupling and not found_halogenation:
            print("Found cross-coupling at late stage but no halogenation in earlier steps")

    return strategy_present
