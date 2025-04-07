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
    Detects a synthetic strategy involving coupling of multiple nitrogen heterocycles
    (pyrazole, pyridone, indole, and other nitrogen heterocycles).
    """
    # Initialize flags
    heterocycles_found = set()
    has_heterocycle_coupling = False

    # Define nitrogen heterocycles to check
    nitrogen_heterocycles = [
        "pyrazole",
        "indole",
        "imidazole",
        "triazole",
        "tetrazole",
        "pyridine",
        "pyrimidine",
        "pyrazine",
        "pyridazine",
        "quinoline",
        "isoquinoline",
        "benzimidazole",
        "indazole",
        "benzotriazole",
    ]

    def dfs_traverse(node, depth=0):
        nonlocal has_heterocycle_coupling, heterocycles_found

        if node["type"] == "mol":
            # Check for heterocycles in molecules
            mol_smiles = node["smiles"]

            # Check for all nitrogen heterocycles
            for heterocycle in nitrogen_heterocycles:
                if checker.check_ring(heterocycle, mol_smiles):
                    heterocycles_found.add(heterocycle)
                    print(f"Found {heterocycle} in molecule: {mol_smiles}")

            # Special check for pyridone (hydroxypyridine or 2-pyridone pattern)
            if checker.check_ring("pyridine", mol_smiles) and (
                checker.check_fg("Enol", mol_smiles)
                or "O=c1ccccn1" in mol_smiles
                or "O=c1ncccc1" in mol_smiles
            ):
                heterocycles_found.add("pyridone")
                print(f"Found pyridone in molecule: {mol_smiles}")

        elif node["type"] == "reaction":
            # Check for heterocycle coupling
            if "metadata" in node and "rsmi" in node["metadata"]:
                rsmi = node["metadata"]["rsmi"]
                reactants = rsmi.split(">")[0].split(".")
                product = rsmi.split(">")[-1]

                # Check if this is a coupling reaction
                coupling_reactions = [
                    "N-arylation (Buchwald-Hartwig/Ullmann-Goldberg)",
                    "Buchwald-Hartwig",
                    "N-arylation_heterocycles",
                    "Ullmann-Goldberg Substitution amine",
                    "Suzuki",
                    "Negishi",
                    "Stille",
                    "Heck_terminal_vinyl",
                    "Sonogashira",
                    "Ullmann condensation",
                    "Goldberg coupling",
                ]

                is_coupling = False
                for reaction_type in coupling_reactions:
                    if checker.check_reaction(reaction_type, rsmi):
                        is_coupling = True
                        print(f"Found coupling reaction {reaction_type}: {rsmi}")
                        break

                # Check for heterocycles in reactants
                reactant_heterocycles = set()
                for reactant in reactants:
                    for heterocycle in nitrogen_heterocycles:
                        if checker.check_ring(heterocycle, reactant):
                            reactant_heterocycles.add(heterocycle)
                            print(f"Found {heterocycle} in reactant: {reactant}")

                    # Special check for pyridone
                    if checker.check_ring("pyridine", reactant) and (
                        checker.check_fg("Enol", reactant)
                        or "O=c1ccccn1" in reactant
                        or "O=c1ncccc1" in reactant
                    ):
                        reactant_heterocycles.add("pyridone")
                        print(f"Found pyridone in reactant: {reactant}")

                # Check for heterocycles in product
                product_heterocycles = set()
                for heterocycle in nitrogen_heterocycles:
                    if checker.check_ring(heterocycle, product):
                        product_heterocycles.add(heterocycle)
                        print(f"Found {heterocycle} in product: {product}")

                # Special check for pyridone in product
                if checker.check_ring("pyridine", product) and (
                    checker.check_fg("Enol", product)
                    or "O=c1ccccn1" in product
                    or "O=c1ncccc1" in product
                ):
                    product_heterocycles.add("pyridone")
                    print(f"Found pyridone in product: {product}")

                print(f"Reactant heterocycles: {reactant_heterocycles}")
                print(f"Product heterocycles: {product_heterocycles}")

                # Check for heterocycle coupling conditions
                if len(reactant_heterocycles) >= 1 and len(product_heterocycles) >= 1:
                    # Case 1: Multiple heterocycles in reactants that are preserved in product
                    if len(
                        reactant_heterocycles
                    ) >= 2 and reactant_heterocycles.issubset(product_heterocycles):
                        has_heterocycle_coupling = True
                        print(
                            f"Detected multiple heterocycle preservation at depth {depth}"
                        )

                    # Case 2: At least one heterocycle in reactants and multiple in product
                    elif len(
                        product_heterocycles
                    ) >= 2 and reactant_heterocycles.intersection(product_heterocycles):
                        has_heterocycle_coupling = True
                        print(f"Detected heterocycle coupling at depth {depth}")

                    # Case 3: Different heterocycles in different reactants
                    elif len(reactants) >= 2:
                        reactant_heterocycle_sets = []
                        for reactant in reactants:
                            r_heterocycles = set()
                            for heterocycle in nitrogen_heterocycles:
                                if checker.check_ring(heterocycle, reactant):
                                    r_heterocycles.add(heterocycle)

                            # Check for pyridone
                            if checker.check_ring("pyridine", reactant) and (
                                checker.check_fg("Enol", reactant)
                                or "O=c1ccccn1" in reactant
                                or "O=c1ncccc1" in reactant
                            ):
                                r_heterocycles.add("pyridone")

                            if r_heterocycles:
                                reactant_heterocycle_sets.append(r_heterocycles)

                        # If we have heterocycles in at least two different reactants
                        if len(reactant_heterocycle_sets) >= 2:
                            # Check if they're different heterocycles
                            all_heterocycles = set()
                            for r_set in reactant_heterocycle_sets:
                                all_heterocycles.update(r_set)

                            if len(all_heterocycles) >= 2:
                                has_heterocycle_coupling = True
                                print(
                                    f"Detected heterocycle coupling from different reactants at depth {depth}"
                                )

                # Even if not a named coupling reaction, check if we're connecting heterocycles
                if not is_coupling and len(reactants) >= 2:
                    # Check if different reactants have different heterocycles
                    # and if the product contains multiple heterocycles
                    if (
                        len(reactant_heterocycles) >= 1
                        and len(product_heterocycles) >= 2
                    ):
                        # Check if product contains all reactant heterocycles
                        if reactant_heterocycles.issubset(product_heterocycles):
                            has_heterocycle_coupling = True
                            print(
                                f"Detected implicit heterocycle coupling at depth {depth}"
                            )

        # Traverse children
        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    # Start traversal from the root
    dfs_traverse(route)

    print(
        f"Final results: heterocycles found = {heterocycles_found}, coupling = {has_heterocycle_coupling}"
    )

    # Return True if we have a heterocycle coupling and at least two different heterocycles
    return has_heterocycle_coupling and len(heterocycles_found) >= 2
