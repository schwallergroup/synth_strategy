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
    Detects a strategy involving SNAr reactions on nitro-activated aromatic rings.
    """
    nitro_activated_snar_found = False

    def dfs_traverse(node):
        nonlocal nitro_activated_snar_found

        if node["type"] == "reaction":
            try:
                # Extract reactants and product
                rsmi = node["metadata"]["rsmi"]
                reactants = rsmi.split(">")[0].split(".")
                product = rsmi.split(">")[-1]

                # Check if this is a nucleophilic aromatic substitution reaction
                if checker.check_reaction(
                    "nucl_sub_aromatic_ortho_nitro", rsmi
                ) or checker.check_reaction("nucl_sub_aromatic_para_nitro", rsmi):
                    print(f"Found potential nitro-activated SNAr reaction: {rsmi}")

                    # Check for nitro-activated haloarene in reactants
                    for reactant in reactants:
                        if checker.check_fg("Nitro group", reactant) and checker.check_fg(
                            "Aromatic halide", reactant
                        ):
                            print(f"Found nitro-activated haloarene: {reactant}")

                            # Check for amine nucleophile in other reactants
                            for other_reactant in reactants:
                                if other_reactant != reactant and (
                                    checker.check_fg("Primary amine", other_reactant)
                                    or checker.check_fg("Secondary amine", other_reactant)
                                    or checker.check_fg("Aniline", other_reactant)
                                ):
                                    print(f"Found amine nucleophile: {other_reactant}")

                                    # Verify product has nitro group and new C-N bond
                                    if checker.check_fg("Nitro group", product) and (
                                        checker.check_fg("Primary amine", product)
                                        or checker.check_fg("Secondary amine", product)
                                        or checker.check_fg("Tertiary amine", product)
                                        or checker.check_fg("Aniline", product)
                                    ):
                                        print(f"Confirmed nitro-activated SNAr product: {product}")
                                        nitro_activated_snar_found = True
                                        return

                # Alternative detection method if reaction type check fails
                if not nitro_activated_snar_found:
                    for reactant in reactants:
                        # Check for nitro-activated haloarene
                        if checker.check_fg("Nitro group", reactant) and checker.check_fg(
                            "Aromatic halide", reactant
                        ):
                            print(f"Found nitro-activated haloarene (alt method): {reactant}")

                            # Check for amine nucleophile
                            amine_found = False
                            for other_reactant in reactants:
                                if other_reactant != reactant and (
                                    checker.check_fg("Primary amine", other_reactant)
                                    or checker.check_fg("Secondary amine", other_reactant)
                                    or checker.check_fg("Aniline", other_reactant)
                                ):
                                    amine_found = True
                                    print(f"Found amine nucleophile (alt method): {other_reactant}")
                                    break

                            # Verify product has nitro group but no halide, and has amine group
                            if amine_found:
                                prod_mol = Chem.MolFromSmiles(product)
                                if (
                                    checker.check_fg("Nitro group", product)
                                    and not checker.check_fg("Aromatic halide", product)
                                    and (
                                        checker.check_fg("Primary amine", product)
                                        or checker.check_fg("Secondary amine", product)
                                        or checker.check_fg("Tertiary amine", product)
                                        or checker.check_fg("Aniline", product)
                                    )
                                ):
                                    print(
                                        f"Confirmed nitro-activated SNAr product (alt method): {product}"
                                    )
                                    nitro_activated_snar_found = True
                                    return
            except Exception as e:
                print(f"Error processing reaction node: {e}")

        # Traverse children
        for child in node.get("children", []):
            dfs_traverse(child)
            if nitro_activated_snar_found:
                return

    # Start traversal
    dfs_traverse(route)

    return nitro_activated_snar_found
