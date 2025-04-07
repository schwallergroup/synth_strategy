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
    Detects if the route contains a hydroboration-oxidation sequence (alkene to alcohol).
    """
    hydroboration_found = False

    def dfs_traverse(node):
        nonlocal hydroboration_found

        if node["type"] == "reaction" and "rsmi" in node.get("metadata", {}):
            rsmi = node["metadata"]["rsmi"]
            reactants = rsmi.split(">")[0].split(".")
            product = rsmi.split(">")[-1]

            # Check for alkene in reactants
            has_alkene = False
            has_oxidant = False

            for reactant in reactants:
                if reactant:
                    mol = Chem.MolFromSmiles(reactant)
                    if mol:
                        alkene_pattern = Chem.MolFromSmarts("[C]=[C]")
                        if mol.HasSubstructMatch(alkene_pattern):
                            has_alkene = True

                        # Check for hydrogen peroxide or similar oxidant
                        if "OO" in reactant:
                            has_oxidant = True

            # Check for alcohol in product
            if has_alkene:
                prod_mol = Chem.MolFromSmiles(product)
                if prod_mol:
                    alcohol_pattern = Chem.MolFromSmarts("[C][OH]")
                    if prod_mol.HasSubstructMatch(alcohol_pattern):
                        print("Found hydroboration-oxidation pattern")
                        hydroboration_found = True

        for child in node.get("children", []):
            dfs_traverse(child)

    dfs_traverse(route)
    return hydroboration_found
