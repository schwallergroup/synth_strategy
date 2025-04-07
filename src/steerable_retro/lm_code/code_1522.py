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
    This function detects if the synthetic route employs multiple different protecting groups
    (e.g., MOM and benzyl) for similar functional groups.
    """
    # Track different protecting groups
    mom_protection = False
    benzyl_protection = False

    def dfs_traverse(node):
        nonlocal mom_protection, benzyl_protection

        if node["type"] == "reaction":
            # Extract reactants and product from reaction SMILES
            rsmi = node["metadata"]["rsmi"]
            reactants_smiles = rsmi.split(">")[0].split(".")
            product_smiles = rsmi.split(">")[-1]

            # Check for MOM group in any molecule
            mom_pattern = Chem.MolFromSmarts("[O]-[CH2]-[O]-[CH3]")

            # Check for benzyl group in any molecule
            benzyl_pattern = Chem.MolFromSmarts("[O]-[CH2]-c1ccccc1")

            # Check reactants and product for MOM
            for smiles in reactants_smiles + [product_smiles]:
                try:
                    mol = Chem.MolFromSmiles(smiles)
                    if mol:
                        if mol.HasSubstructMatch(mom_pattern):
                            mom_protection = True
                        if mol.HasSubstructMatch(benzyl_pattern):
                            benzyl_protection = True
                except:
                    continue

        # Traverse children
        for child in node.get("children", []):
            dfs_traverse(child)

    # Start traversal from the root
    dfs_traverse(route)

    # Return True if both MOM and benzyl protecting groups are found
    return mom_protection and benzyl_protection
