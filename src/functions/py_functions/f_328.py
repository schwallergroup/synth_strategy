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
    This function detects reductive amination (aldehyde + amine -> amine) in the synthetic route.
    """
    reductive_amination_found = False

    def dfs_traverse(node):
        nonlocal reductive_amination_found

        if (
            node["type"] == "reaction"
            and "metadata" in node
            and "rsmi" in node["metadata"]
        ):
            rsmi = node["metadata"]["rsmi"]
            reactants = rsmi.split(">")[0]
            product = rsmi.split(">")[-1]

            # Check if reactants contain aldehyde and amine
            reactants_mol = Chem.MolFromSmiles(reactants)
            if (
                reactants_mol
                and reactants_mol.HasSubstructMatch(Chem.MolFromSmarts("[CH]=O"))
                and reactants_mol.HasSubstructMatch(Chem.MolFromSmarts("[NH]"))
            ):

                # Check if product has new C-N bond
                product_mol = Chem.MolFromSmiles(product)
                if product_mol and product_mol.HasSubstructMatch(
                    Chem.MolFromSmarts("[CH2][N]")
                ):
                    print("Reductive amination detected")
                    reductive_amination_found = True

        for child in node.get("children", []):
            dfs_traverse(child)

    dfs_traverse(route)
    return reductive_amination_found
