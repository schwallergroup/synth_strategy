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
    Detects if the synthesis involves the transformation of a thione (C=S)
    to a thioether (C-S-C) in a heterocyclic system.
    """
    thione_to_thioether = False

    # List of heterocyclic rings to check
    heterocyclic_rings = [
        "thiazole",
        "benzothiazole",
        "thiadiazole",
        "thiophene",
        "benzothiophene",
        "oxathiolane",
        "thiazolidine",
        "isothiazole",
        "imidazole",
        "benzimidazole",
        "pyrimidine",
        "pyrazine",
        "pyridazine",
    ]

    def dfs_traverse(node):
        nonlocal thione_to_thioether

        if node["type"] == "reaction" and "metadata" in node and "rsmi" in node["metadata"]:
            rsmi = node["metadata"]["rsmi"]
            reactants = rsmi.split(">")[0].split(".")
            product = rsmi.split(">")[-1]

            print(f"Analyzing reaction: {rsmi}")

            # Check if this is a thioether formation reaction
            is_thioether_reaction = checker.check_reaction("thioether_nucl_sub", rsmi)
            print(f"Is thioether nucleophilic substitution: {is_thioether_reaction}")

            # Check for thione in reactants
            thione_reactant = None
            for reactant in reactants:
                if checker.check_fg("Thiocarbonyl", reactant):
                    print(f"Found thiocarbonyl in reactant: {reactant}")
                    thione_reactant = reactant
                    break

            # Check for thioether in product
            thioether_in_product = checker.check_fg("Monosulfide", product)
            if thioether_in_product:
                print(f"Found monosulfide in product: {product}")

            # If we have both thione in reactant and thioether in product
            if thione_reactant and thioether_in_product:
                # Check if they are in heterocyclic systems
                reactant_heterocycle = False
                product_heterocycle = False

                for ring in heterocyclic_rings:
                    if checker.check_ring(ring, thione_reactant):
                        print(f"Thiocarbonyl is in a {ring} ring in reactant")
                        reactant_heterocycle = True

                    if checker.check_ring(ring, product):
                        print(f"Monosulfide is in a {ring} ring in product")
                        product_heterocycle = True

                if reactant_heterocycle and product_heterocycle:
                    print("Confirmed thione to thioether transformation in heterocyclic system")
                    thione_to_thioether = True

            # Check for S-alkylation reactions
            if not thione_to_thioether:
                s_alkylation = any(
                    checker.check_reaction(rxn, rsmi)
                    for rxn in [
                        "S-alkylation of thiols",
                        "S-alkylation of thiols (ethyl)",
                        "S-alkylation of thiols with alcohols",
                        "S-alkylation of thiols with alcohols (ethyl)",
                    ]
                )

                if s_alkylation:
                    print(f"Found S-alkylation reaction: {rsmi}")

                    # Check for thione in reactants and thioether in product
                    thione_in_reactants = any(
                        checker.check_fg("Thiocarbonyl", r) for r in reactants
                    )
                    thioether_in_product = checker.check_fg("Monosulfide", product)

                    if thione_in_reactants and thioether_in_product:
                        # Check for heterocyclic systems
                        reactant_heterocycle = any(
                            any(checker.check_ring(ring, r) for ring in heterocyclic_rings)
                            for r in reactants
                        )
                        product_heterocycle = any(
                            checker.check_ring(ring, product) for ring in heterocyclic_rings
                        )

                        if reactant_heterocycle and product_heterocycle:
                            thione_to_thioether = True
                            print(
                                "Confirmed thione to thioether transformation in heterocyclic system via S-alkylation"
                            )

            # Additional check for direct C=S to C-S-C transformation
            if not thione_to_thioether and thione_reactant and thioether_in_product:
                # Check if the reaction involves a thione being converted to a thioether
                # This is a more general check that doesn't rely on specific reaction types
                print("Checking for direct thione to thioether transformation")

                # Check if the sulfur atom from the thione is preserved in the product
                # This is a heuristic check based on the reaction SMILES
                if "=[S:" in rsmi and "[S:" in product and not "=[S:" in product:
                    print("Found potential direct thione to thioether transformation")

                    # Verify it's in a heterocyclic system
                    reactant_heterocycle = any(
                        checker.check_ring(ring, thione_reactant) for ring in heterocyclic_rings
                    )
                    product_heterocycle = any(
                        checker.check_ring(ring, product) for ring in heterocyclic_rings
                    )

                    if reactant_heterocycle and product_heterocycle:
                        thione_to_thioether = True
                        print(
                            "Confirmed direct thione to thioether transformation in heterocyclic system"
                        )

        for child in node.get("children", []):
            dfs_traverse(child)

    dfs_traverse(route)

    print(f"Thione to thioether transformation detected: {thione_to_thioether}")
    return thione_to_thioether
