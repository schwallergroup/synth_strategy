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
    This function detects if the synthetic route employs a benzyl protection followed by
    deprotection sequence for phenols.
    """
    # Track benzyl protection and deprotection events
    benzyl_protection_found = False
    benzyl_deprotection_found = False

    def dfs_traverse(node):
        nonlocal benzyl_protection_found, benzyl_deprotection_found

        if node["type"] == "reaction":
            # Extract reactants and product from reaction SMILES
            rsmi = node["metadata"]["rsmi"]
            reactants_smiles = rsmi.split(">")[0].split(".")
            product_smiles = rsmi.split(">")[-1]

            # Check for benzyl protection: phenol + benzyl source -> benzyl ether
            phenol_pattern = Chem.MolFromSmarts("[OH]-c")
            benzyl_ether_pattern = Chem.MolFromSmarts("[O]-[CH2]-c1ccccc1")

            # Check for protection
            phenol_in_reactants = False
            for reactant in reactants_smiles:
                try:
                    mol = Chem.MolFromSmiles(reactant)
                    if mol and mol.HasSubstructMatch(phenol_pattern):
                        phenol_in_reactants = True
                        break
                except:
                    continue

            if phenol_in_reactants:
                try:
                    product_mol = Chem.MolFromSmiles(product_smiles)
                    if product_mol and product_mol.HasSubstructMatch(benzyl_ether_pattern):
                        print("Found benzyl protection: phenol -> benzyl ether")
                        benzyl_protection_found = True
                except:
                    pass

            # Check for deprotection: benzyl ether -> phenol
            benzyl_in_reactants = False
            for reactant in reactants_smiles:
                try:
                    mol = Chem.MolFromSmiles(reactant)
                    if mol and mol.HasSubstructMatch(benzyl_ether_pattern):
                        benzyl_in_reactants = True
                        break
                except:
                    continue

            if benzyl_in_reactants:
                try:
                    product_mol = Chem.MolFromSmiles(product_smiles)
                    if product_mol and product_mol.HasSubstructMatch(phenol_pattern):
                        print("Found benzyl deprotection: benzyl ether -> phenol")
                        benzyl_deprotection_found = True
                except:
                    pass

        # Traverse children
        for child in node.get("children", []):
            dfs_traverse(child)

    # Start traversal from the root
    dfs_traverse(route)

    # Return True if both protection and deprotection are found
    return benzyl_protection_found and benzyl_deprotection_found
