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
    This function detects aminonitrile formation from a ketone.
    """
    found_aminonitrile_formation = False

    def dfs_traverse(node):
        nonlocal found_aminonitrile_formation

        if node["type"] == "reaction":
            # Extract reactants and product
            rsmi = node["metadata"]["rsmi"]
            reactants_smiles = rsmi.split(">")[0].split(".")
            product_smiles = rsmi.split(">")[-1]

            # Check for aminonitrile formation
            for reactant in reactants_smiles:
                reactant_mol = Chem.MolFromSmiles(reactant)
                product_mol = Chem.MolFromSmiles(product_smiles)

                if reactant_mol is not None and product_mol is not None:
                    # Check for ketone in reactant
                    ketone_pattern = Chem.MolFromSmarts("[C](=[O])[C]")
                    # Check for aminonitrile in product
                    aminonitrile_pattern = Chem.MolFromSmarts("[C]([N])([C]#[N])")

                    if (
                        reactant_mol.HasSubstructMatch(ketone_pattern)
                        and product_mol.HasSubstructMatch(aminonitrile_pattern)
                        and not reactant_mol.HasSubstructMatch(aminonitrile_pattern)
                    ):
                        print("Found aminonitrile formation from ketone")
                        found_aminonitrile_formation = True

        # Continue traversing
        for child in node.get("children", []):
            dfs_traverse(child)

    # Start traversal
    dfs_traverse(route)

    return found_aminonitrile_formation
