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
    This function detects a synthetic strategy involving pyrazole modification
    with a late-stage C-C bond formation involving a nitrile-containing fragment.
    """
    # Track if we found key features
    has_pyrazole = False
    has_bromination = False
    has_n_alkylation = False
    has_late_cc_bond = False
    has_nitrile_fragment = False

    def dfs_traverse(node, depth=0):
        nonlocal has_pyrazole, has_bromination, has_n_alkylation, has_late_cc_bond, has_nitrile_fragment

        if node["type"] == "mol":
            # Check for pyrazole in molecules
            if checker.check_ring("pyrazole", node["smiles"]):
                has_pyrazole = True
                print(f"Found pyrazole in molecule: {node['smiles']}")

            # Check for nitrile groups
            if checker.check_fg("Nitrile", node["smiles"]):
                has_nitrile_fragment = True
                print(f"Found nitrile group in molecule: {node['smiles']}")

            # Check for halomethyl pyrazole which indicates N-alkylation preparation
            if checker.check_ring("pyrazole", node["smiles"]) and (
                "CCl" in node["smiles"]
                or "CBr" in node["smiles"]
                or checker.check_fg("Primary halide", node["smiles"])
            ):
                has_n_alkylation = True
                print(
                    f"Found N-alkylation preparation (halomethyl pyrazole) at depth {depth}: {node['smiles']}"
                )

        elif node["type"] == "reaction" and "metadata" in node and "rsmi" in node["metadata"]:
            rsmi = node["metadata"]["rsmi"]
            reactants = rsmi.split(">")[0].split(".")
            product = rsmi.split(">")[-1]

            # Check for bromination reaction
            if checker.check_reaction("Aromatic bromination", rsmi) or checker.check_reaction(
                "Bromination", rsmi
            ):
                # Verify it's bromination of a pyrazole
                has_pyrazole_reactant = any(checker.check_ring("pyrazole", r) for r in reactants)
                has_pyrazole_product = checker.check_ring("pyrazole", product)
                has_bromine_product = (
                    checker.check_fg("Aromatic halide", product) or "Br" in product
                )

                if has_pyrazole_reactant and has_pyrazole_product and has_bromine_product:
                    has_bromination = True
                    print(f"Found pyrazole bromination reaction at depth {depth}: {rsmi}")

            # Check for N-alkylation
            n_alkylation_reactions = [
                "N-alkylation of primary amines with alkyl halides",
                "N-alkylation of secondary amines with alkyl halides",
                "Alkylation of amines",
            ]

            if any(checker.check_reaction(rxn, rsmi) for rxn in n_alkylation_reactions):
                # Verify it's specifically a pyrazole N-alkylation
                has_pyrazole_reactant = any(checker.check_ring("pyrazole", r) for r in reactants)
                has_pyrazole_product = checker.check_ring("pyrazole", product)

                if has_pyrazole_reactant and has_pyrazole_product:
                    has_n_alkylation = True
                    print(f"Found pyrazole N-alkylation reaction at depth {depth}: {rsmi}")

            # Check for late-stage C-C bond formation with nitrile
            if depth <= 1:  # Allow for depth 0 or 1 for late-stage
                # Common C-C bond forming reactions
                cc_bond_reactions = [
                    "Suzuki",
                    "Negishi",
                    "Stille",
                    "Heck_terminal_vinyl",
                    "Sonogashira",
                    "Kumada cross-coupling",
                    "Hiyama-Denmark Coupling",
                ]

                if any(checker.check_reaction(rxn, rsmi) for rxn in cc_bond_reactions):
                    # Check if one reactant has nitrile and another has pyrazole
                    has_nitrile_reactant = any(checker.check_fg("Nitrile", r) for r in reactants)
                    has_pyrazole_reactant = any(
                        checker.check_ring("pyrazole", r) for r in reactants
                    )

                    if has_nitrile_reactant and has_pyrazole_reactant:
                        has_late_cc_bond = True
                        print(
                            f"Found late-stage C-C bond formation with nitrile at depth {depth}: {rsmi}"
                        )

                # Also check for other C-C bond forming reactions that might not be in our list
                if not has_late_cc_bond:
                    has_nitrile_reactant = any(checker.check_fg("Nitrile", r) for r in reactants)
                    has_pyrazole_reactant = any(
                        checker.check_ring("pyrazole", r) for r in reactants
                    )

                    if (
                        has_nitrile_reactant
                        and has_pyrazole_reactant
                        and checker.check_ring("pyrazole", product)
                    ):
                        # Check if the product has both pyrazole and nitrile
                        if checker.check_fg("Nitrile", product):
                            has_late_cc_bond = True
                            print(
                                f"Found potential late-stage C-C bond formation at depth {depth}: {rsmi}"
                            )

        # Additional check for bromination by examining molecule patterns
        if (
            not has_bromination
            and node["type"] == "mol"
            and checker.check_ring("pyrazole", node["smiles"])
        ):
            if "Br" in node["smiles"] and checker.check_fg("Aromatic halide", node["smiles"]):
                has_bromination = True
                print(f"Found brominated pyrazole molecule at depth {depth}: {node['smiles']}")

        # Traverse children
        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    # Start traversal
    dfs_traverse(route)

    # Check if all required features are present
    strategy_present = (
        has_pyrazole
        and has_bromination
        and has_n_alkylation
        and has_late_cc_bond
        and has_nitrile_fragment
    )

    print(
        f"Pyrazole-based synthesis with late-stage C-C bond strategy detected: {strategy_present}"
    )
    print(f"  - Has pyrazole: {has_pyrazole}")
    print(f"  - Has bromination: {has_bromination}")
    print(f"  - Has N-alkylation: {has_n_alkylation}")
    print(f"  - Has late-stage C-C bond: {has_late_cc_bond}")
    print(f"  - Has nitrile fragment: {has_nitrile_fragment}")

    return strategy_present
