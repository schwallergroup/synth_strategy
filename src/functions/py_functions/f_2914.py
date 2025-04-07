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
    Detects if the synthetic route includes an ester hydrolysis step (ester → carboxylic acid).
    """
    ester_hydrolysis_found = False

    def dfs_traverse(node):
        nonlocal ester_hydrolysis_found

        if (
            node["type"] == "reaction"
            and "metadata" in node
            and "rsmi" in node["metadata"]
        ):
            rsmi = node["metadata"]["rsmi"]
            reactants = rsmi.split(">")[0].split(".")
            product = rsmi.split(">")[-1]

            try:
                reactant_mols = [Chem.MolFromSmiles(r) for r in reactants]
                product_mol = Chem.MolFromSmiles(product)

                # Check if reactant contains ester and product contains carboxylic acid
                has_ester = any(
                    r and r.HasSubstructMatch(Chem.MolFromSmarts("[C](=[O])[O][C]"))
                    for r in reactant_mols
                )
                has_acid = product_mol and product_mol.HasSubstructMatch(
                    Chem.MolFromSmarts("[C](=[O])[OH]")
                )

                if has_ester and has_acid:
                    print("Ester hydrolysis detected")
                    ester_hydrolysis_found = True
            except:
                print("Error processing reaction SMILES for ester hydrolysis detection")

        for child in node.get("children", []):
            dfs_traverse(child)

    dfs_traverse(route)
    return ester_hydrolysis_found
