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

root_data = "/home/andres/Documents/steerable_retro/data"

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


def main(route):
    """
    Detects a synthetic strategy involving aromatic halogenation followed by
    metal-catalyzed coupling.
    """
    halogenation_found = False
    halogenation_depth = -1
    coupling_found = False
    coupling_depth = -1
    halogenated_products = set()  # Track halogenated products

    def dfs_traverse(node, depth=0):
        nonlocal halogenation_found, halogenation_depth, coupling_found, coupling_depth

        if node["type"] == "reaction":
            if "metadata" in node and "rsmi" in node["metadata"]:
                rsmi = node["metadata"]["rsmi"]
                reactants = rsmi.split(">")[0].split(".")
                product = rsmi.split(">")[-1]

                # Check for aromatic halogenation reactions
                halogenation_reactions = [
                    "Aromatic fluorination",
                    "Aromatic chlorination",
                    "Aromatic bromination",
                    "Aromatic iodination",
                ]

                # Check for metal-catalyzed coupling reactions
                coupling_reactions = [
                    "Suzuki coupling with boronic acids",
                    "Suzuki coupling with boronic esters",
                    "Suzuki coupling with boronic acids OTf",
                    "Suzuki coupling with boronic esters OTf",
                    "Negishi coupling",
                    "Stille reaction_aryl",
                    "Stille reaction_vinyl",
                    "Hiyama-Denmark Coupling",
                    "Kumada cross-coupling",
                    "Sonogashira alkyne_aryl halide",
                    "Sonogashira acetylene_aryl halide",
                    "Buchwald-Hartwig",
                    "N-arylation",
                    "Ullmann-Goldberg Substitution amine",
                    "Goldberg coupling",
                    "Heck terminal vinyl",
                    "Heck_terminal_vinyl",
                ]

                # Check if this is a halogenation reaction
                is_halogenation = False
                for rxn_type in halogenation_reactions:
                    if checker.check_reaction(rxn_type, rsmi):
                        is_halogenation = True
                        halogenation_found = True
                        halogenation_depth = depth
                        print(f"Found {rxn_type} at depth {depth}")
                        # Add the product to our tracked halogenated products
                        halogenated_products.add(product)
                        break

                # If not a known halogenation reaction, check for pattern
                if not is_halogenation:
                    # Check if product has aromatic halide but reactants don't
                    product_has_aromatic_halide = checker.check_fg("Aromatic halide", product)
                    reactants_have_aromatic_halide = any(
                        checker.check_fg("Aromatic halide", r) for r in reactants
                    )
                    reactants_have_aromatic = any(
                        checker.check_ring("benzene", r) for r in reactants
                    )

                    if (
                        product_has_aromatic_halide
                        and reactants_have_aromatic
                        and not reactants_have_aromatic_halide
                    ):
                        halogenation_found = True
                        halogenation_depth = depth
                        print(f"Found aromatic halogenation pattern at depth {depth}")
                        halogenated_products.add(product)

                # Check if this is a coupling reaction
                is_coupling = False
                for rxn_type in coupling_reactions:
                    if checker.check_reaction(rxn_type, rsmi):
                        is_coupling = True
                        coupling_found = True
                        coupling_depth = depth
                        print(f"Found {rxn_type} at depth {depth}")
                        break

                # If not a known coupling reaction, check for pattern
                if not is_coupling:
                    # Check if any reactant has aromatic halide
                    reactant_with_aromatic_halide = False
                    for r in reactants:
                        if checker.check_fg("Aromatic halide", r):
                            reactant_with_aromatic_halide = True
                            break

                    # Check if product doesn't have aromatic halide
                    product_has_no_aromatic_halide = not checker.check_fg(
                        "Aromatic halide", product
                    )

                    # Check for coupling pattern: reactant has aromatic halide, product doesn't
                    if reactant_with_aromatic_halide and product_has_no_aromatic_halide:
                        # Additional check: product should have an aromatic ring
                        if checker.check_ring("benzene", product) or any(
                            checker.check_ring(ring, product)
                            for ring in [
                                "pyridine",
                                "pyrimidine",
                                "pyrazine",
                                "pyridazine",
                                "naphthalene",
                                "quinoline",
                                "isoquinoline",
                            ]
                        ):
                            coupling_found = True
                            coupling_depth = depth
                            print(f"Found potential coupling pattern at depth {depth}")

        # Traverse children
        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    # Start traversal
    dfs_traverse(route)

    # Check if the strategy is present
    # In retrosynthesis, higher depth means earlier stage in forward synthesis
    strategy_present = (
        halogenation_found
        and coupling_found
        and halogenation_depth
        > coupling_depth  # In retrosynthesis, halogenation appears after coupling
    )

    if strategy_present:
        print("Found aromatic halogenation followed by coupling strategy")
    else:
        if halogenation_found and coupling_found:
            print(
                f"Found both reactions but in wrong order: halogenation at {halogenation_depth}, coupling at {coupling_depth}"
            )
        elif halogenation_found:
            print("Found only halogenation, no coupling")
        elif coupling_found:
            print("Found only coupling, no halogenation")
        else:
            print("Neither halogenation nor coupling found")

    return strategy_present
