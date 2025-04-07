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
    Detects if the synthesis route involves modifications to heterocyclic cores,
    particularly focusing on quinoline/quinolone transformations.
    """
    has_heterocycle_modification = False

    # List of common heterocyclic rings to check
    heterocycles = [
        "quinoline",
        "isoquinoline",
        "pyridine",
        "pyrrole",
        "furan",
        "thiophene",
        "imidazole",
        "oxazole",
        "thiazole",
        "pyrazole",
        "pyrimidine",
        "pyrazine",
        "indole",
        "benzimidazole",
        "benzoxazole",
        "benzothiazole",
        "purine",
        "piperidine",
        "piperazine",
        "morpholine",
        "thiomorpholine",
    ]

    def dfs_traverse(node, depth=0):
        nonlocal has_heterocycle_modification

        if node["type"] == "reaction":
            # Extract reactants and product
            try:
                rsmi = node["metadata"]["rsmi"]
                reactants_smiles = rsmi.split(">")[0].split(".")
                product_smiles = rsmi.split(">")[-1]

                # Check for heterocycle modifications
                for reactant_smiles in reactants_smiles:
                    # Skip if reactant or product is invalid
                    if not reactant_smiles or not product_smiles:
                        continue

                    # Check for heterocycle transformations
                    reactant_heterocycles = []
                    product_heterocycles = []

                    # Identify heterocycles in reactant
                    for ring in heterocycles:
                        if checker.check_ring(ring, reactant_smiles):
                            reactant_heterocycles.append(ring)
                            print(f"Found {ring} in reactant: {reactant_smiles}")

                    # Identify heterocycles in product
                    for ring in heterocycles:
                        if checker.check_ring(ring, product_smiles):
                            product_heterocycles.append(ring)
                            print(f"Found {ring} in product: {product_smiles}")

                    # Check for specific heterocycle transformations
                    if reactant_heterocycles and product_heterocycles:
                        # Case 1: Different heterocycles in reactant and product
                        if set(reactant_heterocycles) != set(product_heterocycles):
                            has_heterocycle_modification = True
                            print(
                                f"Heterocycle transformation detected: {reactant_heterocycles} → {product_heterocycles}"
                            )

                        # Case 2: Same heterocycle but with functional group modifications
                        else:
                            reactant_mol = Chem.MolFromSmiles(reactant_smiles)
                            product_mol = Chem.MolFromSmiles(product_smiles)

                            if reactant_mol and product_mol:
                                # Check for common functional group modifications on heterocycles
                                fg_groups = [
                                    "Carbonyl",
                                    "Hydroxyl",
                                    "Amine",
                                    "Halide",
                                    "Nitro",
                                ]

                                reactant_fgs = []
                                product_fgs = []

                                for fg in fg_groups:
                                    if checker.check_fg(fg, reactant_smiles):
                                        reactant_fgs.append(fg)
                                    if checker.check_fg(fg, product_smiles):
                                        product_fgs.append(fg)

                                if set(reactant_fgs) != set(product_fgs):
                                    # Check if the reaction involves a known heterocycle modification
                                    if (
                                        checker.check_reaction(
                                            "Friedel-Crafts acylation", rsmi
                                        )
                                        or checker.check_reaction(
                                            "Aromatic nitration", rsmi
                                        )
                                        or checker.check_reaction(
                                            "Aromatic halogenation", rsmi
                                        )
                                        or checker.check_reaction("N-arylation", rsmi)
                                    ):
                                        has_heterocycle_modification = True
                                        print(
                                            f"Functional group modification on heterocycle detected: {rsmi}"
                                        )

                # Specific check for quinoline/quinolone transformations
                if (
                    checker.check_ring("quinoline", reactants_smiles[0])
                    and not checker.check_ring("quinoline", product_smiles)
                    and checker.check_fg("Carbonyl", product_smiles)
                ):
                    has_heterocycle_modification = True
                    print(f"Quinoline to quinolone transformation detected: {rsmi}")

                elif (
                    not checker.check_ring("quinoline", reactants_smiles[0])
                    and checker.check_fg("Carbonyl", reactants_smiles[0])
                    and checker.check_ring("quinoline", product_smiles)
                ):
                    has_heterocycle_modification = True
                    print(f"Quinolone to quinoline transformation detected: {rsmi}")

            except Exception as e:
                print(f"Error processing reaction: {e}")

        # Traverse children
        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    # Start traversal
    dfs_traverse(route)

    print(f"Heterocycle modification strategy: {has_heterocycle_modification}")
    return has_heterocycle_modification
