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
    This function detects if the synthesis route involves a specific sequence of
    functional group interconversions: alkene → ester → alcohol → aldehyde → α-aminonitrile → amide
    """
    # Initialize dictionaries to track transformations and their depths
    transformations = {
        "alkene_to_ester": None,
        "ester_to_alcohol": None,
        "alcohol_to_aldehyde": None,
        "aldehyde_to_aminonitrile": None,
        "aminonitrile_to_amide": None,
    }

    def dfs_traverse(node, depth=0):
        if node["type"] == "reaction":
            if "rsmi" in node.get("metadata", {}):
                rsmi = node["metadata"]["rsmi"]
                reactants = rsmi.split(">")[0].split(".")
                product = rsmi.split(">")[-1]

                # Define patterns for functional groups
                alkene_ester_pattern = Chem.MolFromSmarts("[C]=[C][C](=O)[O][C]")
                ester_pattern = Chem.MolFromSmarts("[C](=O)[O][C]")
                alcohol_pattern = Chem.MolFromSmarts("[CH2][OH]")
                aldehyde_pattern = Chem.MolFromSmarts("[CH]=O")
                aminonitrile_pattern = Chem.MolFromSmarts("[C]([NH2])[C]#[N]")
                amide_pattern = Chem.MolFromSmarts("[NH][C](=O)")

                # Check for alkene to ester transformation
                alkene_in_reactants = any(
                    Chem.MolFromSmiles(r) is not None
                    and Chem.MolFromSmiles(r).HasSubstructMatch(alkene_ester_pattern)
                    for r in reactants
                    if r
                )
                ester_in_product = (
                    Chem.MolFromSmiles(product) is not None
                    and Chem.MolFromSmiles(product).HasSubstructMatch(ester_pattern)
                    and not Chem.MolFromSmiles(product).HasSubstructMatch(
                        alkene_ester_pattern
                    )
                )

                if alkene_in_reactants and ester_in_product:
                    transformations["alkene_to_ester"] = depth
                    print(f"Found alkene to ester transformation at depth {depth}")

                # Check for ester to alcohol transformation
                ester_in_reactants = any(
                    Chem.MolFromSmiles(r) is not None
                    and Chem.MolFromSmiles(r).HasSubstructMatch(ester_pattern)
                    for r in reactants
                    if r
                )
                alcohol_in_product = Chem.MolFromSmiles(
                    product
                ) is not None and Chem.MolFromSmiles(product).HasSubstructMatch(
                    alcohol_pattern
                )

                if ester_in_reactants and alcohol_in_product:
                    transformations["ester_to_alcohol"] = depth
                    print(f"Found ester to alcohol transformation at depth {depth}")

                # Check for alcohol to aldehyde transformation
                alcohol_in_reactants = any(
                    Chem.MolFromSmiles(r) is not None
                    and Chem.MolFromSmiles(r).HasSubstructMatch(alcohol_pattern)
                    for r in reactants
                    if r
                )
                aldehyde_in_product = Chem.MolFromSmiles(
                    product
                ) is not None and Chem.MolFromSmiles(product).HasSubstructMatch(
                    aldehyde_pattern
                )

                if alcohol_in_reactants and aldehyde_in_product:
                    transformations["alcohol_to_aldehyde"] = depth
                    print(f"Found alcohol to aldehyde transformation at depth {depth}")

                # Check for aldehyde to aminonitrile transformation
                aldehyde_in_reactants = any(
                    Chem.MolFromSmiles(r) is not None
                    and Chem.MolFromSmiles(r).HasSubstructMatch(aldehyde_pattern)
                    for r in reactants
                    if r
                )
                aminonitrile_in_product = Chem.MolFromSmiles(
                    product
                ) is not None and Chem.MolFromSmiles(product).HasSubstructMatch(
                    aminonitrile_pattern
                )

                if aldehyde_in_reactants and aminonitrile_in_product:
                    transformations["aldehyde_to_aminonitrile"] = depth
                    print(
                        f"Found aldehyde to aminonitrile transformation at depth {depth}"
                    )

                # Check for aminonitrile to amide transformation
                aminonitrile_in_reactants = any(
                    Chem.MolFromSmiles(r) is not None
                    and Chem.MolFromSmiles(r).HasSubstructMatch(aminonitrile_pattern)
                    for r in reactants
                    if r
                )
                amide_in_product = Chem.MolFromSmiles(
                    product
                ) is not None and Chem.MolFromSmiles(product).HasSubstructMatch(
                    amide_pattern
                )

                if aminonitrile_in_reactants and amide_in_product:
                    transformations["aminonitrile_to_amide"] = depth
                    print(
                        f"Found aminonitrile to amide transformation at depth {depth}"
                    )

        # Continue traversal
        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    # Start traversal from the root
    dfs_traverse(route)

    # Check if we have at least 4 of the 5 transformations in the correct sequence
    valid_transformations = [v for v in transformations.values() if v is not None]
    if len(valid_transformations) >= 4:
        # Check if the transformations are in the correct order
        # In retrosynthetic direction, the depths should be increasing
        # (meaning decreasing in the forward direction)
        ordered_transformations = sorted(
            [(k, v) for k, v in transformations.items() if v is not None],
            key=lambda x: x[1],
        )
        print(f"Transformation sequence: {ordered_transformations}")
        return True

    return False
