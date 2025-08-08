#!/bin/python

"""LM-defined function for strategy description."""

from rdkit.Chem import AllChem, rdFMCS
import copy
from collections import deque
import rdkit.Chem as Chem
from rdkit.Chem import rdMolDescriptors
from rdkit.Chem import rdChemReactions
from rdkit.Chem import AllChem
from rdkit.Chem import rdFMCS
import rdkit.Chem.rdFMCS
from rdkit.Chem import AllChem, Descriptors, rdMolDescriptors
from rdkit import Chem
from rdkit.Chem import Descriptors
from rdkit.Chem import AllChem, rdMolDescriptors
from rdkit.Chem import AllChem, Descriptors, Lipinski
from rdkit.Chem import rdmolops
import re
from rdkit.Chem.Scaffolds import MurckoScaffold
from rdkit.Chem import AllChem, Descriptors
import traceback
import rdkit
from collections import Counter


def main(route):
    """
    This function detects if the synthetic route employs sulfonamide protection
    of a nitrogen atom in an indoline/indole scaffold.
    """
    sulfonamide_protection = False

    def dfs_traverse(node):
        nonlocal sulfonamide_protection

        if node["type"] == "reaction" and "metadata" in node and "rsmi" in node["metadata"]:
            rsmi = node["metadata"]["rsmi"]
            reactants_smiles = rsmi.split(">")[0]
            product_smiles = rsmi.split(">")[-1]

            # Define SMARTS patterns for indoline/indole N and sulfonamide
            indoline_n_pattern = Chem.MolFromSmarts("N1CCc2ccccc12")
            indole_n_pattern = Chem.MolFromSmarts("[nH]1ccc2ccccc12")
            sulfonamide_pattern = Chem.MolFromSmarts("NS(=O)(=O)")

            try:
                reactants_mol = Chem.MolFromSmiles(reactants_smiles)
                product_mol = Chem.MolFromSmiles(product_smiles)

                if reactants_mol and product_mol:
                    # Check if reactants contain indoline/indole N and product contains sulfonamide
                    if (
                        reactants_mol.HasSubstructMatch(indoline_n_pattern)
                        or reactants_mol.HasSubstructMatch(indole_n_pattern)
                    ) and product_mol.HasSubstructMatch(sulfonamide_pattern):
                        print("Detected sulfonamide protection of indoline/indole nitrogen")
                        sulfonamide_protection = True
            except:
                print("Error processing SMILES in sulfonamide protection detection")

        # Continue traversing
        for child in node.get("children", []):
            dfs_traverse(child)

    # Start traversal
    dfs_traverse(route)
    return sulfonamide_protection
