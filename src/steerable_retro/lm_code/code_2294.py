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
    This function detects Knoevenagel-type condensation between aldehyde and active methylene compound.
    """
    knoevenagel_detected = False

    def is_knoevenagel_reaction(rsmi):
        """Helper function to identify Knoevenagel condensation by chemical characteristics"""
        try:
            reactants_part = rsmi.split(">")[0]
            products_part = rsmi.split(">")[-1]

            reactants = reactants_part.split(".")

            # Check if any reactant has an aldehyde or ketone group
            has_carbonyl = False
            has_active_methylene = False

            for reactant in reactants:
                # Check for aldehyde
                if checker.check_fg("Aldehyde", reactant) or checker.check_fg(
                    "Formaldehyde", reactant
                ):
                    has_carbonyl = True
                    print(f"Found aldehyde in reactant: {reactant}")

                # Check for ketone
                if checker.check_fg("Ketone", reactant):
                    has_carbonyl = True
                    print(f"Found ketone in reactant: {reactant}")

                # Check for active methylene compounds (compounds with electron-withdrawing groups)
                # Common active methylene compounds in Knoevenagel: malononitrile, cyanoacetates, malonates
                if "[C]#N" in reactant and "CH2" in reactant:
                    has_active_methylene = True
                    print(f"Found potential active methylene compound: {reactant}")

                # Check for nitrile groups which are common in active methylene compounds
                if checker.check_fg("Nitrile", reactant):
                    # If it has a nitrile and likely a CH2 group adjacent to electron-withdrawing groups
                    if "CH2" in reactant or "[CH2" in reactant:
                        has_active_methylene = True
                        print(f"Found nitrile-containing active methylene compound: {reactant}")

            # Check if product has an α,β-unsaturated structure (C=C-C=O or similar)
            product_mol = Chem.MolFromSmiles(products_part)
            if product_mol and (
                checker.check_fg("Nitrile", products_part) or "C=C" in products_part
            ):
                print(f"Product has characteristics of Knoevenagel product: {products_part}")

                # If we have both required reactant types, this is likely a Knoevenagel condensation
                if has_carbonyl and has_active_methylene:
                    print("Identified as Knoevenagel condensation based on reactants and product")
                    return True

            return False
        except Exception as e:
            print(f"Error in is_knoevenagel_reaction: {e}")
            return False

    def dfs_traverse(node, depth=0):
        nonlocal knoevenagel_detected

        if node["type"] == "reaction":
            if "rsmi" in node.get("metadata", {}):
                rsmi = node["metadata"]["rsmi"]
                print(f"Checking reaction at depth {depth}: {rsmi}")

                # Use the checker function to detect Knoevenagel condensation
                if checker.check_reaction("Knoevenagel Condensation", rsmi):
                    knoevenagel_detected = True
                    print(f"Knoevenagel condensation detected at depth {depth} by checker")

                    # Extract reactants and product for additional verification
                    reactants_smiles = rsmi.split(">")[0].split(".")
                    product_smiles = rsmi.split(">")[-1]

                    # Print details about the detected reaction
                    print(f"Reactants: {reactants_smiles}")
                    print(f"Product: {product_smiles}")
                # Backup method to detect Knoevenagel by chemical characteristics
                elif is_knoevenagel_reaction(rsmi):
                    knoevenagel_detected = True
                    print(
                        f"Knoevenagel condensation detected at depth {depth} by chemical analysis"
                    )

                    # Extract reactants and product for additional verification
                    reactants_smiles = rsmi.split(">")[0].split(".")
                    product_smiles = rsmi.split(">")[-1]

                    # Print details about the detected reaction
                    print(f"Reactants: {reactants_smiles}")
                    print(f"Product: {product_smiles}")

        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    dfs_traverse(route)
    return knoevenagel_detected
