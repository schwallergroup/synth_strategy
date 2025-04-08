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
    Detects if the route contains a Wittig olefination (aldehyde to alkene).
    """
    wittig_found = False

    def dfs_traverse(node):
        nonlocal wittig_found

        if node["type"] == "reaction" and "rsmi" in node.get("metadata", {}):
            rsmi = node["metadata"]["rsmi"]
            reactants = rsmi.split(">")[0].split(".")
            product = rsmi.split(">")[-1]

            # Check for aldehyde in reactants
            has_aldehyde = False
            has_phosphonium = False

            for reactant in reactants:
                if reactant:
                    mol = Chem.MolFromSmiles(reactant)
                    if mol:
                        aldehyde_pattern = Chem.MolFromSmarts("[C]=O")
                        if mol.HasSubstructMatch(aldehyde_pattern):
                            has_aldehyde = True

                        # Check for phosphonium salt or ylide (simplified)
                        if "P+" in reactant or "P" in reactant:
                            has_phosphonium = True

            # Check for alkene in product
            if has_aldehyde and has_phosphonium:
                prod_mol = Chem.MolFromSmiles(product)
                if prod_mol:
                    alkene_pattern = Chem.MolFromSmarts("[C]=[C]")
                    if prod_mol.HasSubstructMatch(alkene_pattern):
                        print("Found Wittig olefination")
                        wittig_found = True

        for child in node.get("children", []):
            dfs_traverse(child)

    dfs_traverse(route)
    return wittig_found
