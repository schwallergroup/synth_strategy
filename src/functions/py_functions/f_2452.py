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
    Detects Cbz (benzyloxycarbonyl) protection of an amine.
    """
    found_cbz_protection = False

    def dfs_traverse(node):
        nonlocal found_cbz_protection

        if node["type"] == "reaction":
            if "metadata" in node and "rsmi" in node["metadata"]:
                rsmi = node["metadata"]["rsmi"]
                reactants = rsmi.split(">")[0].split(".")
                product = rsmi.split(">")[-1]

                # Check for amine to Cbz-protected amine conversion
                amine_pattern = Chem.MolFromSmarts("[NH2]")
                cbz_pattern = Chem.MolFromSmarts("[NH][C](=[O])[O][CH2][c]")

                # Check if reactants contain chloroformate and amine
                has_amine = False
                has_chloroformate = False

                for reactant in reactants:
                    reactant_mol = Chem.MolFromSmiles(reactant)
                    if reactant_mol:
                        if reactant_mol.HasSubstructMatch(amine_pattern):
                            has_amine = True
                        if reactant_mol.HasSubstructMatch(
                            Chem.MolFromSmarts("[Cl][C](=[O])[O][CH2][c]")
                        ):
                            has_chloroformate = True

                # Check if product contains Cbz group
                product_mol = Chem.MolFromSmiles(product)
                if (
                    product_mol
                    and product_mol.HasSubstructMatch(cbz_pattern)
                    and has_amine
                    and has_chloroformate
                ):
                    found_cbz_protection = True
                    print("Found Cbz protection of amine")

        # Traverse children
        for child in node.get("children", []):
            dfs_traverse(child)

    # Start traversal
    dfs_traverse(route)

    return found_cbz_protection
