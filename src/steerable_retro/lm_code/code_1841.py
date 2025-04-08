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
    Detects if the synthesis route includes esterification of an alcohol.
    """
    esterification_found = False

    def dfs_traverse(node):
        nonlocal esterification_found

        if node["type"] == "reaction" and "metadata" in node and "rsmi" in node["metadata"]:
            rsmi = node["metadata"]["rsmi"]
            reactants = rsmi.split(">")[0].split(".")
            product = rsmi.split(">")[-1]

            # Check for alcohol to ester conversion
            for reactant in reactants:
                reactant_mol = Chem.MolFromSmiles(reactant)
                product_mol = Chem.MolFromSmiles(product)

                if reactant_mol and product_mol:
                    alcohol_patt = Chem.MolFromSmarts("[C][OH]")
                    ester_patt = Chem.MolFromSmarts("[C][O][C](=O)[C]")

                    if (
                        reactant_mol.HasSubstructMatch(alcohol_patt)
                        and product_mol.HasSubstructMatch(ester_patt)
                        and not reactant_mol.HasSubstructMatch(ester_patt)
                    ):
                        # Check if another reactant is an acylating agent
                        for r in reactants:
                            if "C(=O)O" in r or "C(=O)Cl" in r:
                                print(f"Alcohol esterification detected: {rsmi}")
                                esterification_found = True
                                break

        for child in node.get("children", []):
            dfs_traverse(child)

    dfs_traverse(route)
    return esterification_found
