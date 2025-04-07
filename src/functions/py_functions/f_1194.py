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


def main(route):
    """
    This function detects if the synthetic route involves sequential modification
    of a side chain (alkene → alcohol → leaving group → nitrile) on an aromatic core.
    """
    # Track if we've seen each transformation
    alkene_to_alcohol = False
    alcohol_to_leaving_group = False
    leaving_group_to_nitrile = False

    def dfs_traverse(node, depth=0):
        nonlocal alkene_to_alcohol, alcohol_to_leaving_group, leaving_group_to_nitrile

        if node["type"] == "reaction":
            if "rsmi" in node.get("metadata", {}):
                rsmi = node["metadata"]["rsmi"]
                reactants = rsmi.split(">")[0].split(".")
                product = rsmi.split(">")[-1]

                # Create patterns for functional groups
                alkene_pattern = Chem.MolFromSmarts("C=C")
                alcohol_pattern = Chem.MolFromSmarts("[C][OH]")
                mesylate_pattern = Chem.MolFromSmarts("[O][S](=[O])(=[O])[C]")
                nitrile_pattern = Chem.MolFromSmarts("C#N")

                product_mol = Chem.MolFromSmiles(product)

                # Check for transformations by comparing reactants and products
                for reactant in reactants:
                    reactant_mol = Chem.MolFromSmiles(reactant)

                    if not reactant_mol or not product_mol:
                        continue

                    # Check for alkene to alcohol transformation
                    if (
                        reactant_mol.HasSubstructMatch(alkene_pattern)
                        and product_mol.HasSubstructMatch(alcohol_pattern)
                        and not reactant_mol.HasSubstructMatch(alcohol_pattern)
                    ):
                        alkene_to_alcohol = True
                        print(
                            f"Alkene to alcohol transformation found at depth {depth}"
                        )

                    # Check for alcohol to leaving group transformation
                    if (
                        reactant_mol.HasSubstructMatch(alcohol_pattern)
                        and product_mol.HasSubstructMatch(mesylate_pattern)
                        and not reactant_mol.HasSubstructMatch(mesylate_pattern)
                    ):
                        alcohol_to_leaving_group = True
                        print(
                            f"Alcohol to leaving group transformation found at depth {depth}"
                        )

                    # Check for leaving group to nitrile transformation
                    if (
                        reactant_mol.HasSubstructMatch(mesylate_pattern)
                        and product_mol.HasSubstructMatch(nitrile_pattern)
                        and not reactant_mol.HasSubstructMatch(nitrile_pattern)
                    ):
                        leaving_group_to_nitrile = True
                        print(
                            f"Leaving group to nitrile transformation found at depth {depth}"
                        )

        # Continue traversal
        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    # Start traversal
    dfs_traverse(route)

    # Return True if we found at least two of the three transformations in sequence
    transformations_found = sum(
        [alkene_to_alcohol, alcohol_to_leaving_group, leaving_group_to_nitrile]
    )
    return transformations_found >= 2
