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
    Detects if the synthesis route involves acylation of an aromatic amine.
    """
    acylation_found = False

    def dfs_traverse(node):
        nonlocal acylation_found

        if node["type"] == "reaction" and node.get("metadata", {}).get("rsmi"):
            try:
                rsmi = node["metadata"]["rsmi"]
                print(f"Examining reaction: {rsmi}")

                # Check for acylation reaction types directly
                acylation_reactions = [
                    "Acylation of Nitrogen Nucleophiles by Acyl/Thioacyl/Carbamoyl Halides and Analogs_N",
                    "Schotten-Baumann to ester",
                    "Acyl chloride with primary amine to amide (Schotten-Baumann)",
                    "Acylation of primary amines",
                    "Acylation of secondary amines",
                    "Schotten-Baumann_amide",
                ]

                for reaction_type in acylation_reactions:
                    if checker.check_reaction(reaction_type, rsmi):
                        print(f"Found reaction type: {reaction_type}")

                        # Correctly split the reaction SMILES
                        parts = rsmi.split(">")
                        if len(parts) >= 3:
                            reactants_smiles = parts[0]
                            products_smiles = parts[2]

                            # Check if reactants contain aromatic amine
                            reactants = reactants_smiles.split(".")
                            for reactant in reactants:
                                # Check for aniline or primary amine on benzene
                                if checker.check_fg("Aniline", reactant) or (
                                    checker.check_fg("Primary amine", reactant)
                                    and checker.check_ring("benzene", reactant)
                                ):
                                    print(
                                        f"Found aromatic amine in reactant: {reactant}"
                                    )

                                    # Check if product contains acylated amine
                                    products = products_smiles.split(".")
                                    for product in products:
                                        # Check for amide connected to aromatic ring
                                        if (
                                            checker.check_fg("Primary amide", product)
                                            or checker.check_fg(
                                                "Secondary amide", product
                                            )
                                            or checker.check_fg(
                                                "Tertiary amide", product
                                            )
                                        ) and checker.check_ring("benzene", product):
                                            print(
                                                f"Found acylated aromatic amine in product: {product}"
                                            )
                                            acylation_found = True
                                            return

                # Additional check for specific cases where the reaction might not be directly classified
                parts = rsmi.split(">")
                if len(parts) >= 3:
                    reactants_smiles = parts[0]
                    products_smiles = parts[2]

                    reactants = reactants_smiles.split(".")
                    products = products_smiles.split(".")

                    # Check for aromatic amine in reactants and amide in products
                    aromatic_amine_in_reactants = False
                    for reactant in reactants:
                        if checker.check_fg("Aniline", reactant) or (
                            checker.check_fg("Primary amine", reactant)
                            and checker.check_ring("benzene", reactant)
                        ):
                            aromatic_amine_in_reactants = True
                            print(
                                f"Found aromatic amine in reactant (additional check): {reactant}"
                            )
                            break

                    if aromatic_amine_in_reactants:
                        for product in products:
                            if (
                                (
                                    checker.check_fg("Primary amide", product)
                                    or checker.check_fg("Secondary amide", product)
                                    or checker.check_fg("Tertiary amide", product)
                                )
                                and checker.check_ring("benzene", product)
                                and not checker.check_fg("Aniline", product)
                            ):
                                print(
                                    f"Found acylated aromatic amine in product (additional check): {product}"
                                )
                                acylation_found = True
                                return

            except Exception as e:
                print(f"Error processing reaction: {e}")

        # Continue traversal
        for child in node.get("children", []):
            dfs_traverse(child)

    # Start traversal from the root
    dfs_traverse(route)
    print(f"Acylation found: {acylation_found}")
    return acylation_found
