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
    This function detects if the synthesis involves late-stage isoxazole formation
    as the final step in a multi-ring system synthesis.
    """
    isoxazole_formed = False
    depth_of_formation = None

    def count_rings(mol_smiles):
        """Count the number of rings in a molecule"""
        mol = Chem.MolFromSmiles(mol_smiles)
        if not mol:
            return 0
        return mol.GetRingInfo().NumRings()

    def contains_isoxazole_pattern(smiles):
        """Check if the molecule contains an isoxazole pattern"""
        # This is a backup check in case the checker.check_ring function doesn't catch it
        mol = Chem.MolFromSmiles(smiles)
        if not mol:
            return False
        # Look for the isoxazole pattern: 5-membered ring with O and two N atoms
        isoxazole_pattern = Chem.MolFromSmarts("c1noc[n,c]1")
        return mol.HasSubstructMatch(isoxazole_pattern)

    def dfs_traverse(node, depth=0):
        nonlocal isoxazole_formed, depth_of_formation

        if node["type"] == "reaction" and "metadata" in node and "rsmi" in node["metadata"]:
            rsmi = node["metadata"]["rsmi"]
            reactants = rsmi.split(">")[0].split(".")
            product = rsmi.split(">")[-1]

            print(f"Analyzing reaction at depth {depth}: {rsmi}")

            # Check if product contains isoxazole
            has_isoxazole_in_product = checker.check_ring("isoxazole", product)
            backup_isoxazole_check = contains_isoxazole_pattern(product)

            print(
                f"Product isoxazole check: checker={has_isoxazole_in_product}, backup={backup_isoxazole_check}"
            )

            if has_isoxazole_in_product or backup_isoxazole_check:
                print(f"Product contains isoxazole: {product}")

                # Check if any reactant contains isoxazole
                reactants_have_isoxazole = any(
                    checker.check_ring("isoxazole", reactant)
                    or contains_isoxazole_pattern(reactant)
                    for reactant in reactants
                )

                print(f"Reactants have isoxazole: {reactants_have_isoxazole}")

                if not reactants_have_isoxazole:
                    print(f"Isoxazole formation detected at depth {depth}")

                    # Check if this is a multi-ring system
                    product_ring_count = count_rings(product)
                    print(f"Product has {product_ring_count} rings")

                    # Check for various isoxazole formation reactions
                    huisgen_1_3 = checker.check_reaction("Huisgen 1,3 dipolar cycloaddition", rsmi)
                    huisgen_alkene = checker.check_reaction(
                        "Huisgen alkene-azide 1,3 dipolar cycloaddition", rsmi
                    )
                    huisgen_alkyne = checker.check_reaction(
                        "Huisgen alkyne-azide 1,3 dipolar cycloaddition", rsmi
                    )
                    huisgen_cu = checker.check_reaction("{Huisgen_Cu-catalyzed_1,4-subst}", rsmi)
                    huisgen_ru = checker.check_reaction("{Huisgen_Ru-catalyzed_1,5_subst}", rsmi)
                    huisgen_disubst = checker.check_reaction("{Huisgen_disubst-alkyne}", rsmi)

                    print(
                        f"Reaction checks: Huisgen_1_3={huisgen_1_3}, alkene={huisgen_alkene}, alkyne={huisgen_alkyne}, Cu={huisgen_cu}, Ru={huisgen_ru}, disubst={huisgen_disubst}"
                    )

                    is_isoxazole_formation = (
                        huisgen_1_3
                        or huisgen_alkene
                        or huisgen_alkyne
                        or huisgen_cu
                        or huisgen_ru
                        or huisgen_disubst
                    )

                    # If none of the specific reactions match, check if it's any kind of isoxazole formation
                    if not is_isoxazole_formation:
                        # If we detected isoxazole in product but not in reactants, and it's not one of the known reactions,
                        # we'll assume it's still an isoxazole formation reaction
                        print(
                            "No specific isoxazole formation reaction detected, but isoxazole appears to be formed"
                        )
                        is_isoxazole_formation = True

                    if product_ring_count > 1 and is_isoxazole_formation:
                        isoxazole_formed = True
                        depth_of_formation = depth
                        print(f"Multi-ring isoxazole formation confirmed at depth {depth}")

        # Traverse children
        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    # Start traversal
    dfs_traverse(route)

    # Check if isoxazole formation is late-stage (depth = 1)
    if isoxazole_formed and depth_of_formation == 1:
        print("Late-stage isoxazole formation in multi-ring system confirmed")
        return True

    # Also check for depth 0 as it might be the final step in some cases
    if isoxazole_formed and depth_of_formation == 0:
        print("Late-stage isoxazole formation in multi-ring system confirmed (at target molecule)")
        return True

    print(f"Late-stage isoxazole formation not detected. Formation depth: {depth_of_formation}")
    return False
