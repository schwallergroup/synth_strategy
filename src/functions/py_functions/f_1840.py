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
    Detects if the synthesis route includes oxidation of an alcohol to a ketone.
    """
    oxidation_found = False

    def dfs_traverse(node):
        nonlocal oxidation_found

        if (
            node["type"] == "reaction"
            and "metadata" in node
            and "rsmi" in node["metadata"]
        ):
            rsmi = node["metadata"]["rsmi"]
            reactants = rsmi.split(">")[0].split(".")
            product = rsmi.split(">")[-1]

            # Check for alcohol to ketone conversion
            for reactant in reactants:
                reactant_mol = Chem.MolFromSmiles(reactant)
                product_mol = Chem.MolFromSmiles(product)

                if reactant_mol and product_mol:
                    alcohol_patt = Chem.MolFromSmarts("[C][OH]")
                    ketone_patt = Chem.MolFromSmarts("[C](=O)[C]")

                    if (
                        reactant_mol.HasSubstructMatch(alcohol_patt)
                        and product_mol.HasSubstructMatch(ketone_patt)
                        and not reactant_mol.HasSubstructMatch(ketone_patt)
                    ):
                        print(f"Alcohol oxidation detected: {rsmi}")
                        oxidation_found = True

        for child in node.get("children", []):
            dfs_traverse(child)

    dfs_traverse(route)
    return oxidation_found
