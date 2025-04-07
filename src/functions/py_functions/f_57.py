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
    This function detects if the route contains at least 3 distinct functional group interconversions.
    """
    transformations = set()

    def dfs_traverse(node):
        nonlocal transformations

        if (
            node["type"] == "reaction"
            and "metadata" in node
            and "rsmi" in node["metadata"]
        ):
            rsmi = node["metadata"]["rsmi"]
            reactants = rsmi.split(">")[0].split(".")
            product = rsmi.split(">")[-1]

            # Define patterns for functional groups
            patterns = {
                "nitrile": "[#6]C#N",
                "amine": "[#6][NH2]",
                "alcohol": "[#6][OH]",
                "ester": "[#6]C(=O)O[#6]",
                "aldehyde": "[#6][CH]=O",
                "carboxylic_acid": "[#6]C(=O)[OH]",
                "amide": "[#6]C(=O)N[#6]",
                "phthalimide": "[#6]N1C(=O)c2ccccc2C1=O",
            }

            # Check reactants and products for functional groups
            reactant_groups = set()
            for reactant in reactants:
                reactant_mol = Chem.MolFromSmiles(reactant)
                if reactant_mol:
                    for group_name, pattern in patterns.items():
                        if reactant_mol.HasSubstructMatch(Chem.MolFromSmarts(pattern)):
                            reactant_groups.add(group_name)

            product_mol = Chem.MolFromSmiles(product)
            product_groups = set()
            if product_mol:
                for group_name, pattern in patterns.items():
                    if product_mol.HasSubstructMatch(Chem.MolFromSmarts(pattern)):
                        product_groups.add(group_name)

            # Identify transformations
            for r_group in reactant_groups:
                for p_group in product_groups:
                    if r_group != p_group:
                        transformation = f"{r_group}_to_{p_group}"
                        transformations.add(transformation)
                        print(f"Detected transformation: {transformation}")

        # Traverse children
        for child in node.get("children", []):
            dfs_traverse(child)

    # Start traversal
    dfs_traverse(route)

    # Return True if at least 3 distinct transformations were found
    return len(transformations) >= 3
