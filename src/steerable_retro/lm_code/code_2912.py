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
    Detects if the synthetic route includes a nitro reduction step (NO2 → NH2).
    """
    nitro_reduction_found = False

    def dfs_traverse(node):
        nonlocal nitro_reduction_found

        if node["type"] == "reaction" and "metadata" in node and "rsmi" in node["metadata"]:
            rsmi = node["metadata"]["rsmi"]
            reactants = rsmi.split(">")[0].split(".")
            product = rsmi.split(">")[-1]

            # Check if any reactant contains nitro group and product contains amine at same position
            try:
                reactant_mols = [Chem.MolFromSmiles(r) for r in reactants]
                product_mol = Chem.MolFromSmiles(product)

                if any(
                    reactant_mol.HasSubstructMatch(Chem.MolFromSmarts("[N+](=[O])[O-]"))
                    for reactant_mol in reactant_mols
                    if reactant_mol
                ):
                    if product_mol and product_mol.HasSubstructMatch(Chem.MolFromSmarts("[NH2]")):
                        print("Nitro reduction detected")
                        nitro_reduction_found = True
            except:
                print("Error processing reaction SMILES for nitro reduction detection")

        for child in node.get("children", []):
            dfs_traverse(child)

    dfs_traverse(route)
    return nitro_reduction_found
