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
    This function detects the use of malonate chemistry for carbon framework construction.
    """
    malonate_chemistry_used = False

    def dfs_traverse(node):
        nonlocal malonate_chemistry_used

        if (
            node["type"] == "reaction"
            and "metadata" in node
            and "rsmi" in node["metadata"]
        ):
            rsmi = node["metadata"]["rsmi"]
            reactants = rsmi.split(">")[0].split(".")
            product = rsmi.split(">")[-1]

            # Check for malonate pattern in reactants
            for reactant in reactants:
                try:
                    reactant_mol = Chem.MolFromSmiles(reactant)
                    # SMARTS pattern for malonate or similar diester
                    malonate_pattern = Chem.MolFromSmarts("C(=O)OC.C(=O)OC")

                    if reactant_mol and reactant_mol.HasSubstructMatch(
                        malonate_pattern
                    ):
                        # Check if this is used for C-C bond formation
                        product_mol = Chem.MolFromSmiles(product)
                        if product_mol:
                            print(
                                "Detected malonate-based carbon framework construction"
                            )
                            malonate_chemistry_used = True
                            break
                except:
                    pass

        for child in node.get("children", []):
            dfs_traverse(child)

    dfs_traverse(route)
    return malonate_chemistry_used
