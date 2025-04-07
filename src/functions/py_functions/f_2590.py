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
    This function detects a strategy involving the assembly of a complex molecule with multiple heteroaromatic rings.
    """
    # Track if we've found each ring system
    has_isoxazole = False
    has_benzothiophene = False
    has_benzodioxole = False

    # Track if we've found reactions that form these rings
    isoxazole_formation = False
    benzothiophene_formation = False
    benzodioxole_formation = False

    # Check if the final product contains all three rings
    final_product_smiles = route["smiles"]

    # Check for isoxazole
    final_has_isoxazole = checker.check_ring("isoxazole", final_product_smiles)
    print(f"Final product has isoxazole: {final_has_isoxazole}")

    # Check for benzothiophene - either directly or as a benzene fused to thiophene
    final_has_benzothiophene = checker.check_ring(
        "benzothiophene", final_product_smiles
    ) or (
        checker.check_ring("benzene", final_product_smiles)
        and checker.check_ring("thiophene", final_product_smiles)
    )
    print(f"Final product has benzothiophene: {final_has_benzothiophene}")

    # Check for benzodioxole structure - specifically looking for dioxolane fused to benzene
    final_has_benzodioxole = checker.check_ring(
        "dioxolane", final_product_smiles
    ) and checker.check_ring("benzene", final_product_smiles)

    # If dioxolane check fails, look for the OCO pattern in the SMILES which indicates 1,3-benzodioxole
    if not final_has_benzodioxole and "OCO" in final_product_smiles:
        final_has_benzodioxole = True

    print(f"Final product has benzodioxole-like structure: {final_has_benzodioxole}")

    final_product_has_all_rings = (
        final_has_isoxazole and final_has_benzothiophene and final_has_benzodioxole
    )

    if final_product_has_all_rings:
        print(f"Final product contains all three ring systems: {final_product_smiles}")

    def dfs_traverse(node, depth=0):
        nonlocal has_isoxazole, has_benzothiophene, has_benzodioxole
        nonlocal isoxazole_formation, benzothiophene_formation, benzodioxole_formation

        if node["type"] == "mol":
            mol_smiles = node["smiles"]

            # Check for isoxazole
            if checker.check_ring("isoxazole", mol_smiles):
                has_isoxazole = True
                print(f"Found isoxazole at depth {depth}: {mol_smiles}")

            # Check for benzothiophene - either directly or as a benzene fused to thiophene
            if checker.check_ring("benzothiophene", mol_smiles) or (
                checker.check_ring("benzene", mol_smiles)
                and checker.check_ring("thiophene", mol_smiles)
            ):
                has_benzothiophene = True
                print(f"Found benzothiophene at depth {depth}: {mol_smiles}")

            # Check for benzodioxole structure - specifically looking for dioxolane fused to benzene
            if checker.check_ring("dioxolane", mol_smiles) and checker.check_ring(
                "benzene", mol_smiles
            ):
                has_benzodioxole = True
                print(
                    f"Found benzodioxole-like structure at depth {depth}: {mol_smiles}"
                )
            # If dioxolane check fails, look for the OCO pattern in the SMILES which indicates 1,3-benzodioxole
            elif "OCO" in mol_smiles and checker.check_ring("benzene", mol_smiles):
                has_benzodioxole = True
                print(
                    f"Found benzodioxole-like structure at depth {depth}: {mol_smiles}"
                )

        elif (
            node["type"] == "reaction"
            and "metadata" in node
            and "rsmi" in node["metadata"]
        ):
            rxn_smiles = node["metadata"]["rsmi"]
            reactants = rxn_smiles.split(">")[0].split(".")
            product = rxn_smiles.split(">")[-1]

            # Check for reactions that form isoxazole
            if checker.check_reaction("Formation of NOS Heterocycles", rxn_smiles):
                # Check if this reaction produces an isoxazole that wasn't in the reactants
                if checker.check_ring("isoxazole", product) and not any(
                    checker.check_ring("isoxazole", r) for r in reactants
                ):
                    isoxazole_formation = True
                    print(
                        f"Found isoxazole formation reaction at depth {depth}: {rxn_smiles}"
                    )

            # Check for reactions that form benzothiophene
            if checker.check_reaction(
                "benzothiophene", rxn_smiles
            ) or checker.check_reaction("thiophene", rxn_smiles):
                # Check if this reaction produces a benzothiophene that wasn't in the reactants
                has_benzothiophene_product = checker.check_ring(
                    "benzothiophene", product
                ) or (
                    checker.check_ring("benzene", product)
                    and checker.check_ring("thiophene", product)
                )
                has_benzothiophene_reactants = any(
                    checker.check_ring("benzothiophene", r)
                    or (
                        checker.check_ring("benzene", r)
                        and checker.check_ring("thiophene", r)
                    )
                    for r in reactants
                )

                if has_benzothiophene_product and not has_benzothiophene_reactants:
                    benzothiophene_formation = True
                    print(
                        f"Found benzothiophene formation reaction at depth {depth}: {rxn_smiles}"
                    )

            # Check for reactions that form benzodioxole
            if checker.check_reaction(
                "Williamson Ether Synthesis", rxn_smiles
            ) or checker.check_reaction("{Williamson ether}", rxn_smiles):
                # Check if this reaction produces a benzodioxole-like structure that wasn't in the reactants
                has_benzodioxole_product = (
                    checker.check_ring("dioxolane", product)
                    and checker.check_ring("benzene", product)
                ) or ("OCO" in product and checker.check_ring("benzene", product))
                has_benzodioxole_reactants = any(
                    (
                        checker.check_ring("dioxolane", r)
                        and checker.check_ring("benzene", r)
                    )
                    or ("OCO" in r and checker.check_ring("benzene", r))
                    for r in reactants
                )

                if has_benzodioxole_product and not has_benzodioxole_reactants:
                    benzodioxole_formation = True
                    print(
                        f"Found benzodioxole formation reaction at depth {depth}: {rxn_smiles}"
                    )

        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    dfs_traverse(route)

    # We need to find all three ring systems and the final product should have all rings
    ring_systems_found = has_isoxazole and has_benzothiophene and has_benzodioxole
    ring_formation_reactions = (
        isoxazole_formation or benzothiophene_formation or benzodioxole_formation
    )

    print(f"Ring systems found: {ring_systems_found}")
    print(f"Ring formation reactions found: {ring_formation_reactions}")
    print(f"Final product has all rings: {final_product_has_all_rings}")

    # If we found all ring systems and the final product has all rings, consider it a success
    # Even if we didn't detect specific ring formation reactions
    return ring_systems_found and final_product_has_all_rings
