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
    Detects if the synthesis involves heterocycles throughout the synthesis.
    Focuses on common heterocycles like benzothiazole, benzoxazole, benzimidazole, etc.
    """
    # List of heterocycles to check
    heterocycles = [
        "benzothiazole",
        "benzoxazole",
        "benzimidazole",
        "indole",
        "quinoline",
        "isoquinoline",
        "pyridine",
        "pyrimidine",
        "pyrazine",
        "pyridazine",
        "furan",
        "thiophene",
        "pyrrole",
        "imidazole",
        "oxazole",
        "thiazole",
        "triazole",
        "tetrazole",
        "isoxazole",
        "isothiazole",
        "oxadiazole",
        "thiadiazole",
    ]

    has_heterocycle = False
    heterocycle_found = None

    def dfs_traverse(node, depth=0):
        nonlocal has_heterocycle, heterocycle_found

        if node["type"] == "mol":
            # Check for heterocycles in molecule
            for heterocycle in heterocycles:
                if checker.check_ring(heterocycle, node["smiles"]):
                    has_heterocycle = True
                    heterocycle_found = heterocycle
                    print(
                        f"Found {heterocycle} heterocycle in molecule at depth {depth}: {node['smiles']}"
                    )
                    break

        elif node["type"] == "reaction" and "metadata" in node and "rsmi" in node["metadata"]:
            # Check if this is a heterocycle formation reaction
            rsmi = node["metadata"]["rsmi"]
            product = rsmi.split(">")[-1]
            reactants = rsmi.split(">")[0].split(".")

            # Check if product contains any heterocycle
            product_heterocycles = []
            for heterocycle in heterocycles:
                if checker.check_ring(heterocycle, product):
                    product_heterocycles.append(heterocycle)
                    has_heterocycle = True
                    heterocycle_found = heterocycle
                    print(f"Found {heterocycle} in product at depth {depth}: {product}")

            # For each heterocycle found in product, check if it's formed in this reaction
            for heterocycle in product_heterocycles:
                # Check if any reactant contains this heterocycle
                reactant_has_heterocycle = any(
                    checker.check_ring(heterocycle, r) for r in reactants if r
                )

                # If no reactant has this heterocycle but product does, it's a formation reaction
                if not reactant_has_heterocycle:
                    print(f"Found {heterocycle} formation reaction at depth {depth}: {rsmi}")

                    # Check specific formation reactions based on the heterocycle
                    if heterocycle == "benzothiazole":
                        if checker.check_reaction("benzothiazole formation from aldehyde", rsmi):
                            print(f"Identified as benzothiazole formation from aldehyde")
                        elif checker.check_reaction(
                            "benzothiazole formation from acyl halide", rsmi
                        ):
                            print(f"Identified as benzothiazole formation from acyl halide")
                        elif checker.check_reaction(
                            "benzothiazole formation from ester/carboxylic acid", rsmi
                        ):
                            print(
                                f"Identified as benzothiazole formation from ester/carboxylic acid"
                            )
                        elif checker.check_reaction("{benzothiazole}", rsmi):
                            print(f"Identified as benzothiazole formation (generic)")

                    elif heterocycle == "benzoxazole":
                        if checker.check_reaction("benzoxazole formation from aldehyde", rsmi):
                            print(f"Identified as benzoxazole formation from aldehyde")
                        elif checker.check_reaction("benzoxazole formation from acyl halide", rsmi):
                            print(f"Identified as benzoxazole formation from acyl halide")
                        elif checker.check_reaction(
                            "benzoxazole formation from ester/carboxylic acid", rsmi
                        ):
                            print(f"Identified as benzoxazole formation from ester/carboxylic acid")
                        elif checker.check_reaction("benzoxazole formation (intramolecular)", rsmi):
                            print(f"Identified as benzoxazole formation (intramolecular)")
                        elif checker.check_reaction("{benzoxazole_arom-aldehyde}", rsmi):
                            print(f"Identified as benzoxazole formation from aromatic aldehyde")
                        elif checker.check_reaction("{benzoxazole_carboxylic-acid}", rsmi):
                            print(f"Identified as benzoxazole formation from carboxylic acid")

                    elif heterocycle == "benzimidazole":
                        if checker.check_reaction("benzimidazole formation from aldehyde", rsmi):
                            print(f"Identified as benzimidazole formation from aldehyde")
                        elif checker.check_reaction(
                            "benzimidazole formation from acyl halide", rsmi
                        ):
                            print(f"Identified as benzimidazole formation from acyl halide")
                        elif checker.check_reaction(
                            "benzimidazole formation from ester/carboxylic acid", rsmi
                        ):
                            print(
                                f"Identified as benzimidazole formation from ester/carboxylic acid"
                            )
                        elif checker.check_reaction(
                            "{benzimidazole_derivatives_carboxylic-acid/ester}", rsmi
                        ):
                            print(
                                f"Identified as benzimidazole formation from carboxylic acid/ester"
                            )
                        elif checker.check_reaction("{benzimidazole_derivatives_aldehyde}", rsmi):
                            print(f"Identified as benzimidazole formation from aldehyde")

        # Traverse children
        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    # Start traversal from root
    dfs_traverse(route)

    if has_heterocycle:
        print(f"Heterocycle-containing strategy detected: {heterocycle_found}")
    else:
        print("No heterocycle-containing strategy detected")

    return has_heterocycle
