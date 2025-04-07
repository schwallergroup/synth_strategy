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
    This function detects if the synthetic route involves esterification.
    """
    esterification_found = False

    def dfs_traverse(node):
        nonlocal esterification_found

        if (
            node["type"] == "reaction"
            and "metadata" in node
            and "rsmi" in node["metadata"]
        ):
            rsmi = node["metadata"]["rsmi"]
            reactants = rsmi.split(">")[0].split(".")
            product = rsmi.split(">")[-1]

            # Check if reactants contain a carboxylic acid and an alcohol, and product contains an ester
            reactant_mols = [Chem.MolFromSmiles(r) for r in reactants if r]
            product_mol = Chem.MolFromSmiles(product) if product else None

            if product_mol and len(reactant_mols) >= 2:
                # Check for carboxylic acid in reactants
                acid_pattern = Chem.MolFromSmarts("[CX3](=[OX1])[OX2H]")
                has_acid = any(
                    mol.HasSubstructMatch(acid_pattern) for mol in reactant_mols if mol
                )

                # Check for alcohol in reactants
                alcohol_pattern = Chem.MolFromSmarts("[OX2H]")
                has_alcohol = any(
                    mol.HasSubstructMatch(alcohol_pattern)
                    for mol in reactant_mols
                    if mol
                )

                # Check for ester in product
                ester_pattern = Chem.MolFromSmarts("[CX3](=[OX1])[OX2][C]")
                has_ester = (
                    product_mol.HasSubstructMatch(ester_pattern)
                    if product_mol
                    else False
                )

                if has_acid and has_alcohol and has_ester:
                    print("Esterification detected")
                    esterification_found = True

        for child in node.get("children", []):
            dfs_traverse(child)

    dfs_traverse(route)
    return esterification_found
