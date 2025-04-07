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
    This function detects multiple amine functional group interconversions in the synthesis.
    It tracks transformations between nitro, amine, amide, and isothiocyanate groups.
    """
    amine_pattern = Chem.MolFromSmarts("c[NH2]")
    nitro_pattern = Chem.MolFromSmarts("c[N+](=O)[O-]")
    amide_pattern = Chem.MolFromSmarts("c[NH]C(=O)")
    isothiocyanate_pattern = Chem.MolFromSmarts("c[N]=C=S")

    amine_transformations = 0

    def dfs_traverse(node):
        nonlocal amine_transformations

        if node["type"] == "reaction":
            # Extract reactants and product
            rsmi = node["metadata"]["rsmi"]
            reactants_smiles = rsmi.split(">")[0].split(".")
            product_smiles = rsmi.split(">")[-1]

            product_mol = Chem.MolFromSmiles(product_smiles)
            if not product_mol:
                return

            # Check for various transformations
            for reactant_smiles in reactants_smiles:
                reactant_mol = Chem.MolFromSmiles(reactant_smiles)
                if not reactant_mol:
                    continue

                # Nitro to amine
                if (
                    reactant_mol.HasSubstructMatch(nitro_pattern)
                    and product_mol.HasSubstructMatch(amine_pattern)
                    and not reactant_mol.HasSubstructMatch(amine_pattern)
                ):
                    print("Nitro to amine transformation detected")
                    amine_transformations += 1

                # Amine to amide
                if (
                    reactant_mol.HasSubstructMatch(amine_pattern)
                    and product_mol.HasSubstructMatch(amide_pattern)
                    and not reactant_mol.HasSubstructMatch(amide_pattern)
                ):
                    print("Amine to amide transformation detected")
                    amine_transformations += 1

                # Isothiocyanate to amine
                if (
                    reactant_mol.HasSubstructMatch(isothiocyanate_pattern)
                    and product_mol.HasSubstructMatch(amine_pattern)
                    and not reactant_mol.HasSubstructMatch(amine_pattern)
                ):
                    print("Isothiocyanate to amine transformation detected")
                    amine_transformations += 1

        # Traverse children
        for child in node.get("children", []):
            dfs_traverse(child)

    # Start traversal from the root
    dfs_traverse(route)
    return amine_transformations >= 2  # At least 2 different amine transformations
