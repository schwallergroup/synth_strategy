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
    This function detects if the synthesis follows a linear peptide synthesis strategy
    with at least 3 amino acid incorporations.
    """
    amino_acid_incorporations = 0

    # SMARTS patterns
    amino_acid_pattern = Chem.MolFromSmarts("[NH2][CH]([*])[C](=[O])[OH]")
    amide_pattern = Chem.MolFromSmarts("[C](=[O])[NH][CH]([*])[C]")

    def dfs_traverse(node):
        nonlocal amino_acid_incorporations

        if node["type"] == "reaction":
            if "metadata" in node and "rsmi" in node["metadata"]:
                rsmi = node["metadata"]["rsmi"]
                reactants = rsmi.split(">")[0].split(".")
                product = rsmi.split(">")[-1]

                # Check if one of the reactants is an amino acid
                for reactant in reactants:
                    mol = Chem.MolFromSmiles(reactant)
                    if mol and mol.HasSubstructMatch(amino_acid_pattern):
                        # Check if the product has a new amide bond
                        product_mol = Chem.MolFromSmiles(product)
                        if product_mol and product_mol.HasSubstructMatch(amide_pattern):
                            amino_acid_incorporations += 1
                            print(
                                f"Amino acid incorporation detected in reaction: {rsmi}"
                            )
                            break

        # Process children
        for child in node.get("children", []):
            dfs_traverse(child)

    # Start traversal
    dfs_traverse(route)

    result = (
        amino_acid_incorporations >= 2
    )  # At least 2 amino acid incorporations for a tripeptide
    print(f"Linear peptide synthesis strategy: {result}")
    print(f"Amino acid incorporations: {amino_acid_incorporations}")

    return result
