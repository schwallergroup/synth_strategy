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
    Detects if the synthesis route contains a silyl protection/deprotection sequence
    for a primary alcohol.
    """
    protection_found = False
    deprotection_found = False

    def dfs_traverse(node, depth=0):
        nonlocal protection_found, deprotection_found

        if node["type"] == "reaction":
            if "rsmi" in node.get("metadata", {}):
                rsmi = node["metadata"]["rsmi"]
                reactants = rsmi.split(">")[0].split(".")
                product = rsmi.split(">")[-1]

                # Check for silyl protection (alcohol + silyl chloride -> silyl ether)
                if not protection_found:
                    has_alcohol = False
                    has_silyl_chloride = False

                    for reactant in reactants:
                        mol = Chem.MolFromSmiles(reactant)
                        if mol:
                            if mol.HasSubstructMatch(Chem.MolFromSmarts("[CX4][OX2H]")):
                                has_alcohol = True
                            if mol.HasSubstructMatch(Chem.MolFromSmarts("[Cl][Si]")):
                                has_silyl_chloride = True

                    prod_mol = Chem.MolFromSmiles(product)
                    if prod_mol and prod_mol.HasSubstructMatch(
                        Chem.MolFromSmarts("[CX4][OX2][Si]")
                    ):
                        if has_alcohol and has_silyl_chloride:
                            protection_found = True
                            print(f"Found silyl protection at depth {depth}")

                # Check for silyl deprotection (silyl ether -> alcohol)
                if not deprotection_found:
                    for reactant in reactants:
                        mol = Chem.MolFromSmiles(reactant)
                        if mol and mol.HasSubstructMatch(
                            Chem.MolFromSmarts("[CX4][OX2][Si]")
                        ):
                            prod_mol = Chem.MolFromSmiles(product)
                            if prod_mol and prod_mol.HasSubstructMatch(
                                Chem.MolFromSmarts("[CX4][OX2H]")
                            ):
                                deprotection_found = True
                                print(f"Found silyl deprotection at depth {depth}")

        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    dfs_traverse(route)
    return protection_found and deprotection_found
