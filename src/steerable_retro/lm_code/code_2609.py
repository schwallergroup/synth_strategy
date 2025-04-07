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
    Detects if the synthesis involves the formation of a diaryl ether linkage.
    """
    diaryl_ether_formed = False

    def dfs_traverse(node):
        nonlocal diaryl_ether_formed

        if node["type"] == "reaction" and "metadata" in node and "rsmi" in node["metadata"]:
            rsmi = node["metadata"]["rsmi"]
            reactants = rsmi.split(">")[0].split(".")
            product = rsmi.split(">")[-1]

            # Check for phenol in reactants
            phenol_found = False
            for reactant in reactants:
                react_mol = Chem.MolFromSmiles(reactant)
                if react_mol:
                    phenol_patt = Chem.MolFromSmarts("[c]-[OH]")
                    if react_mol.HasSubstructMatch(phenol_patt):
                        phenol_found = True
                        break

            # Check for diaryl ether in product
            prod_mol = Chem.MolFromSmiles(product)
            if prod_mol:
                diaryl_ether_patt = Chem.MolFromSmarts("[c]-[O]-[c]")
                if prod_mol.HasSubstructMatch(diaryl_ether_patt) and phenol_found:
                    diaryl_ether_formed = True

        for child in node.get("children", []):
            dfs_traverse(child)

    dfs_traverse(route)

    print(f"Diaryl ether formation detected: {diaryl_ether_formed}")
    return diaryl_ether_formed
