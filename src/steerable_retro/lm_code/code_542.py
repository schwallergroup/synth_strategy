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
    Detects if the synthesis involves reduction of an amide to an amine.
    """
    amide_reduction_detected = False

    def dfs_traverse(node):
        nonlocal amide_reduction_detected

        if node["type"] == "reaction":
            # Extract reactants and product
            rsmi = node["metadata"]["rsmi"]
            reactants = rsmi.split(">")[0].split(".")
            product = rsmi.split(">")[-1]

            try:
                prod_mol = Chem.MolFromSmiles(product)
                react_mols = [Chem.MolFromSmiles(r) for r in reactants if r]

                # Check for amide in reactants
                amide_patt = Chem.MolFromSmarts("[C;!$(C=C)](=O)[N]")

                # Check if any reactant has an amide
                amide_in_reactants = any(
                    m and m.HasSubstructMatch(amide_patt) for m in react_mols if m
                )

                # Check if product has a CH2-N where the amide was
                if amide_in_reactants and prod_mol:
                    # This is a simplified check - in a real implementation, you'd need to track
                    # the specific atoms involved in the transformation
                    amine_patt = Chem.MolFromSmarts("[CH2][N;!$(NC=O)]")
                    if prod_mol.HasSubstructMatch(amine_patt):
                        print("Amide reduction detected")
                        amide_reduction_detected = True
            except:
                print("Error in processing reaction SMILES")

        # Process children
        for child in node.get("children", []):
            dfs_traverse(child)

    # Start traversal
    dfs_traverse(route)

    return amide_reduction_detected
