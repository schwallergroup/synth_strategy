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
    This function detects if the route includes a bromination step of a heterocycle,
    particularly a thiophene.
    """
    found = False

    # List of heterocycles to check
    heterocycles = [
        "thiophene",
        "furan",
        "pyrrole",
        "imidazole",
        "oxazole",
        "thiazole",
        "pyridine",
        "pyrimidine",
        "pyrazine",
        "pyridazine",
        "triazole",
        "tetrazole",
    ]

    def dfs_traverse(node):
        nonlocal found

        if node["type"] == "reaction" and "rsmi" in node.get("metadata", {}):
            rsmi = node["metadata"]["rsmi"]
            reactants = rsmi.split(">")[0].split(".")
            product = rsmi.split(">")[-1]

            # Check if this is a bromination reaction
            is_bromination = (
                checker.check_reaction("Aromatic bromination", rsmi)
                or checker.check_reaction("Bromination", rsmi)
                or checker.check_reaction(
                    "Wohl-Ziegler bromination benzyl primary", rsmi
                )
                or checker.check_reaction(
                    "Wohl-Ziegler bromination benzyl secondary", rsmi
                )
                or checker.check_reaction(
                    "Wohl-Ziegler bromination benzyl tertiary", rsmi
                )
                or checker.check_reaction(
                    "Wohl-Ziegler bromination allyl primary", rsmi
                )
                or checker.check_reaction(
                    "Wohl-Ziegler bromination allyl secondary", rsmi
                )
                or checker.check_reaction(
                    "Wohl-Ziegler bromination allyl tertiary", rsmi
                )
            )

            # If it's a bromination reaction, check for heterocycles
            if is_bromination:
                print(f"Found bromination reaction: {rsmi}")

                # Check for heterocycle in reactants
                reactant_with_heterocycle = None
                heterocycle_found = None

                for reactant in reactants:
                    if not reactant:
                        continue

                    for heterocycle in heterocycles:
                        if checker.check_ring(heterocycle, reactant):
                            heterocycle_found = heterocycle
                            reactant_with_heterocycle = reactant
                            print(f"Found {heterocycle} in reactant: {reactant}")
                            break

                    if heterocycle_found:
                        break

                # If a heterocycle was found in reactants, check if it was brominated
                if heterocycle_found and reactant_with_heterocycle:
                    # Check if product has the same heterocycle
                    if checker.check_ring(heterocycle_found, product):
                        # Check if product has aromatic halide with bromine
                        if (
                            checker.check_fg("Aromatic halide", product)
                            and "Br" in product
                        ):
                            # Count bromines in heterocycle reactant and product
                            reactant_br_count = reactant_with_heterocycle.count("Br")
                            product_br_count = product.count("Br")

                            # If product has more bromines than the heterocycle reactant
                            if product_br_count > reactant_br_count:
                                found = True
                                print(
                                    f"Found bromination of {heterocycle_found}: {rsmi}"
                                )
                                # Give special attention to thiophene as mentioned in function description
                                if heterocycle_found == "thiophene":
                                    print(
                                        f"Specifically found bromination of thiophene!"
                                    )

            # If not identified as bromination by reaction type, check for bromination patterns
            elif not is_bromination:
                # Look for heterocycle in reactants
                reactant_with_heterocycle = None
                heterocycle_found = None

                for reactant in reactants:
                    if not reactant:
                        continue

                    for heterocycle in heterocycles:
                        if checker.check_ring(heterocycle, reactant):
                            heterocycle_found = heterocycle
                            reactant_with_heterocycle = reactant
                            print(f"Found {heterocycle} in reactant: {reactant}")
                            break

                    if heterocycle_found:
                        break

                # If heterocycle found, check for bromination pattern
                if heterocycle_found and reactant_with_heterocycle:
                    # Check if any reactant contains bromine (potential brominating agent)
                    has_bromine_reagent = any(
                        "Br" in r for r in reactants if r != reactant_with_heterocycle
                    )

                    # Check if product has the same heterocycle
                    if (
                        checker.check_ring(heterocycle_found, product)
                        and has_bromine_reagent
                    ):
                        # Check if product has aromatic halide with bromine
                        if (
                            checker.check_fg("Aromatic halide", product)
                            and "Br" in product
                        ):
                            # Count bromines in heterocycle reactant and product
                            reactant_br_count = reactant_with_heterocycle.count("Br")
                            product_br_count = product.count("Br")

                            # If product has more bromines than the heterocycle reactant
                            if product_br_count > reactant_br_count:
                                found = True
                                print(
                                    f"Found bromination of {heterocycle_found} by pattern matching: {rsmi}"
                                )
                                # Give special attention to thiophene as mentioned in function description
                                if heterocycle_found == "thiophene":
                                    print(
                                        f"Specifically found bromination of thiophene!"
                                    )

        # Traverse children
        for child in node.get("children", []):
            dfs_traverse(child)

    # Start traversal
    dfs_traverse(route)

    print(f"Bromination of heterocycle: {found}")
    return found
