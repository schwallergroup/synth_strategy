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
    Detects if the synthesis route uses TMS protection/deprotection
    sequence for an alkyne.
    """
    tms_protection = False
    tms_deprotection = False

    def dfs_traverse(node):
        nonlocal tms_protection, tms_deprotection

        if node["type"] == "reaction":
            if "rsmi" in node.get("metadata", {}):
                rsmi = node["metadata"]["rsmi"]
                reactants = rsmi.split(">")[0].split(".")
                product = rsmi.split(">")[-1]

                # Check for TMS protection (terminal alkyne to TMS-alkyne)
                terminal_alkyne_pattern = Chem.MolFromSmarts("[C]#[CH]")
                tms_alkyne_pattern = Chem.MolFromSmarts("[C]#[C][Si]([C])([C])[C]")

                # Check for TMS deprotection (TMS-alkyne to terminal alkyne)
                if any(
                    Chem.MolFromSmiles(r)
                    and Chem.MolFromSmiles(r).HasSubstructMatch(terminal_alkyne_pattern)
                    for r in reactants
                    if Chem.MolFromSmiles(r)
                ):
                    prod_mol = Chem.MolFromSmiles(product)
                    if prod_mol and prod_mol.HasSubstructMatch(tms_alkyne_pattern):
                        print("TMS protection detected")
                        tms_protection = True

                if any(
                    Chem.MolFromSmiles(r)
                    and Chem.MolFromSmiles(r).HasSubstructMatch(tms_alkyne_pattern)
                    for r in reactants
                    if Chem.MolFromSmiles(r)
                ):
                    prod_mol = Chem.MolFromSmiles(product)
                    if prod_mol and prod_mol.HasSubstructMatch(terminal_alkyne_pattern):
                        print("TMS deprotection detected")
                        tms_deprotection = True

        # Traverse children
        for child in node.get("children", []):
            dfs_traverse(child)

    dfs_traverse(route)
    return (
        tms_protection and tms_deprotection
    )  # Both protection and deprotection must be detected
