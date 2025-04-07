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
    Detects formation of an ester with a tertiary alcohol containing two thiophene rings.
    """
    found_bis_thiophene_ester = False

    def dfs_traverse(node):
        nonlocal found_bis_thiophene_ester

        if node["type"] == "reaction" and "metadata" in node and "rsmi" in node["metadata"]:
            rsmi = node["metadata"]["rsmi"]
            reactants_part = rsmi.split(">")[0]
            products_part = rsmi.split(">")[-1]

            print(f"Examining reaction: {rsmi}")

            # Check for esterification or related reactions
            is_esterification = checker.check_reaction("Esterification of Carboxylic Acids", rsmi)
            is_hydrolysis = checker.check_reaction(
                "Hydrolysis or Hydrogenolysis of Carboxylic Esters or Thioesters", rsmi
            )

            if is_esterification or is_hydrolysis:
                print(f"Found {'esterification' if is_esterification else 'hydrolysis'} reaction")

                # For esterification: alcohol in reactants, ester in products
                # For hydrolysis: ester in reactants, alcohol in products
                alcohol_molecules = (
                    reactants_part.split(".") if is_esterification else products_part.split(".")
                )
                ester_molecules = (
                    products_part.split(".") if is_esterification else reactants_part.split(".")
                )

                # Look for tertiary alcohol with two thiophene rings
                for mol_smiles in alcohol_molecules:
                    print(f"Checking molecule for tertiary alcohol: {mol_smiles}")

                    if checker.check_fg("Tertiary alcohol", mol_smiles):
                        print("Found tertiary alcohol")

                        # Check for thiophene rings
                        if checker.check_ring("thiophene", mol_smiles):
                            thiophene_indices = checker.get_ring_atom_indices(
                                "thiophene", mol_smiles
                            )
                            thiophene_count = len(thiophene_indices)
                            print(f"Found {thiophene_count} thiophene rings")

                            if thiophene_count >= 2:
                                print("Found tertiary alcohol with at least 2 thiophene rings")

                                # Check if there's a corresponding ester
                                for ester_smiles in ester_molecules:
                                    if checker.check_fg("Ester", ester_smiles):
                                        print(
                                            "Found ester in the corresponding part of the reaction"
                                        )
                                        found_bis_thiophene_ester = True
                                        return  # Early return once found

            # Also check for other reactions that might form this specific ester
            for reactant in reactants_part.split("."):
                if checker.check_fg("Tertiary alcohol", reactant) and checker.check_ring(
                    "thiophene", reactant
                ):
                    thiophene_indices = checker.get_ring_atom_indices("thiophene", reactant)
                    if len(thiophene_indices) >= 2:
                        print("Found tertiary alcohol with at least 2 thiophene rings in reactants")

                        for product in products_part.split("."):
                            if checker.check_fg("Ester", product) and checker.check_ring(
                                "thiophene", product
                            ):
                                product_thiophene_indices = checker.get_ring_atom_indices(
                                    "thiophene", product
                                )
                                if len(product_thiophene_indices) >= 2:
                                    print("Found ester with at least 2 thiophene rings in products")
                                    found_bis_thiophene_ester = True
                                    return  # Early return once found

        # Continue traversing if not found yet
        for child in node.get("children", []):
            dfs_traverse(child)

    # Start traversal
    dfs_traverse(route)
    return found_bis_thiophene_ester
