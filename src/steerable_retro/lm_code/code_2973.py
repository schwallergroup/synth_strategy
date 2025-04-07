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
    Detects if the route incorporates a heterocyclic component (like pyrimidine)
    in a late-stage coupling reaction.
    """
    has_heterocycle_incorporation = False

    # List of heterocycles to check
    heterocycles = [
        "pyrimidine",
        "pyridine",
        "pyrazole",
        "imidazole",
        "oxazole",
        "thiazole",
        "furan",
        "thiophene",
        "pyrrole",
        "isoxazole",
        "isothiazole",
        "triazole",
        "tetrazole",
        "oxadiazole",
        "thiadiazole",
        "benzoxazole",
        "benzothiazole",
        "benzimidazole",
        "indole",
        "quinoline",
        "isoquinoline",
    ]

    # List of coupling reactions to check
    coupling_reactions = [
        "Suzuki coupling with boronic acids",
        "Suzuki coupling with boronic esters",
        "Suzuki coupling with boronic acids OTf",
        "Suzuki coupling with boronic esters OTf",
        "Suzuki coupling with sulfonic esters",
        "Stille reaction_aryl",
        "Negishi coupling",
        "Heck terminal vinyl",
        "Buchwald-Hartwig/Ullmann-Goldberg/N-arylation primary amine",
        "Buchwald-Hartwig/Ullmann-Goldberg/N-arylation secondary amine",
        "Sonogashira acetylene_aryl halide",
        "Sonogashira alkyne_aryl halide",
        "N-arylation (Buchwald-Hartwig/Ullmann-Goldberg)",
        "Ullmann condensation",
        "Hiyama-Denmark Coupling",
        "Kumada cross-coupling",
        "Aryllithium cross-coupling",
        "Suzuki",
        "{Suzuki}",
        "Stille",
        "{Stille}",
        "Negishi",
        "{Negishi}",
        "Buchwald-Hartwig",
        "{Buchwald-Hartwig}",
    ]

    def dfs_traverse(node, depth=0):
        nonlocal has_heterocycle_incorporation

        print(f"Checking node at depth {depth}: {node['type']}")

        if (
            node["type"] == "reaction" and depth <= 3
        ):  # Late stage (including slightly earlier reactions)
            try:
                rsmi = node["metadata"].get("rsmi", "")
                if not rsmi:
                    print(f"No reaction SMILES at depth {depth}")
                    return

                print(f"Examining reaction at depth {depth}: {rsmi}")
                reactants = rsmi.split(">")[0].split(".")
                product = rsmi.split(">")[-1]

                # Check if it's a coupling reaction
                is_coupling = False
                for reaction_type in coupling_reactions:
                    print(f"Checking if reaction is {reaction_type}")
                    if checker.check_reaction(reaction_type, rsmi):
                        is_coupling = True
                        print(f"Found coupling reaction: {reaction_type} at depth {depth}")
                        break

                # Manual check for Suzuki coupling if the standard check failed
                if (
                    not is_coupling
                    and "OB(" in rsmi
                    and ("Cl" in rsmi or "Br" in rsmi or "I" in rsmi)
                ):
                    print("Manual detection of Suzuki coupling with boronic ester")
                    is_coupling = True

                if not is_coupling:
                    print(f"Not a coupling reaction at depth {depth}")
                    return

                # Check for heterocycles in reactants
                heterocycle_reactant = None
                other_reactants = []
                heterocycle_in_reactant = None

                for reactant in reactants:
                    has_heterocycle = False
                    for heterocycle in heterocycles:
                        if checker.check_ring(heterocycle, reactant):
                            has_heterocycle = True
                            heterocycle_in_reactant = heterocycle
                            print(f"Found {heterocycle} in reactant: {reactant}")
                            break

                    if has_heterocycle:
                        if heterocycle_reactant is None:
                            heterocycle_reactant = reactant
                        else:
                            other_reactants.append(reactant)
                    else:
                        other_reactants.append(reactant)

                # If no heterocycle found in reactants, return
                if heterocycle_reactant is None:
                    print(f"No heterocycle found in reactants at depth {depth}")
                    return

                # Check if product has the heterocycle
                product_has_heterocycle = False
                detected_heterocycle = None
                for heterocycle in heterocycles:
                    if checker.check_ring(heterocycle, product):
                        product_has_heterocycle = True
                        detected_heterocycle = heterocycle
                        print(f"Found {heterocycle} in product: {product}")
                        break

                if product_has_heterocycle:
                    # Check if the heterocycle was actually incorporated (not just modified)
                    main_scaffold_has_heterocycle = False
                    for reactant in other_reactants:
                        if detected_heterocycle and checker.check_ring(
                            detected_heterocycle, reactant
                        ):
                            main_scaffold_has_heterocycle = True
                            print(
                                f"Heterocycle {detected_heterocycle} already present in other reactant"
                            )
                            break

                    if not main_scaffold_has_heterocycle:
                        # If the heterocycle in product matches the one in the heterocycle_reactant
                        if detected_heterocycle == heterocycle_in_reactant:
                            has_heterocycle_incorporation = True
                            print(f"Confirmed heterocycle incorporation at depth {depth}")
                        else:
                            print(
                                f"Heterocycle in product ({detected_heterocycle}) doesn't match reactant ({heterocycle_in_reactant})"
                            )
                else:
                    print(f"Product does not contain a heterocycle at depth {depth}")

            except Exception as e:
                print(f"Error processing reaction node: {e}")

        # Continue traversal
        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    # Start traversal
    dfs_traverse(route)
    return has_heterocycle_incorporation
