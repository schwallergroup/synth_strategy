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
    Detects if the synthetic route includes an esterification step (carboxylic acid → ester).
    """
    esterification_found = False

    def dfs_traverse(node):
        nonlocal esterification_found

        if node["type"] == "reaction" and "metadata" in node and "rsmi" in node["metadata"]:
            rsmi = node["metadata"]["rsmi"]
            reactants = rsmi.split(">")[0].split(".")
            product = rsmi.split(">")[-1]

            try:
                reactant_mols = [Chem.MolFromSmiles(r) for r in reactants]
                product_mol = Chem.MolFromSmiles(product)

                # Check if reactants contain carboxylic acid and alcohol, and product contains ester
                has_acid = any(
                    r and r.HasSubstructMatch(Chem.MolFromSmarts("[C](=[O])[OH]"))
                    for r in reactant_mols
                )
                has_alcohol = any(
                    r and r.HasSubstructMatch(Chem.MolFromSmarts("[OH][C]")) for r in reactant_mols
                )
                has_ester = product_mol and product_mol.HasSubstructMatch(
                    Chem.MolFromSmarts("[C](=[O])[O][C]")
                )

                if has_acid and has_alcohol and has_ester:
                    print("Esterification detected")
                    esterification_found = True
            except:
                print("Error processing reaction SMILES for esterification detection")

        for child in node.get("children", []):
            dfs_traverse(child)

    dfs_traverse(route)
    return esterification_found
