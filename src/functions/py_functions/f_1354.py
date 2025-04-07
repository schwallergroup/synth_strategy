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
    Detects a strategy involving cyclization of a terminal alkene to form a ring.
    """
    found_cyclization = False

    def dfs_traverse(node):
        nonlocal found_cyclization

        if node["type"] == "reaction":
            rsmi = node["metadata"]["rsmi"]
            reactants_smiles = rsmi.split(">")[0]
            product_smiles = rsmi.split(">")[-1]

            # Check for terminal alkene in reactants
            try:
                reactant = Chem.MolFromSmiles(reactants_smiles)
                product = Chem.MolFromSmiles(product_smiles)

                if reactant and product:
                    # Check for terminal alkene
                    terminal_alkene_pattern = Chem.MolFromSmarts("C=C[#1,C]")
                    has_terminal_alkene = reactant.HasSubstructMatch(
                        terminal_alkene_pattern
                    )

                    # Count rings in reactants and product
                    reactant_rings = len(Chem.GetSSSR(reactant))
                    product_rings = len(Chem.GetSSSR(product))

                    if has_terminal_alkene and product_rings > reactant_rings:
                        print("Found terminal alkene cyclization")
                        found_cyclization = True
            except:
                pass

        # Traverse children
        for child in node.get("children", []):
            dfs_traverse(child)

    # Start traversal
    dfs_traverse(route)

    return found_cyclization
