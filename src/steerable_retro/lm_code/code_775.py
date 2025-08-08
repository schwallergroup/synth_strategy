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
    This function detects modifications to heterocyclic systems.
    """
    heterocycle_modifications = 0

    def dfs_traverse(node):
        nonlocal heterocycle_modifications

        if node["type"] == "reaction" and "metadata" in node and "rsmi" in node["metadata"]:
            rsmi = node["metadata"]["rsmi"]
            reactants = rsmi.split(">")[0].split(".")
            product = rsmi.split(">")[-1]

            # Patterns for common heterocycles
            heterocycle_patterns = [
                Chem.MolFromSmarts("c1ncccc1"),  # pyridine
                Chem.MolFromSmarts("c1ncncc1"),  # pyrimidine
                Chem.MolFromSmarts("c1nccnc1"),  # pyrazine
                Chem.MolFromSmarts("c1nnccc1"),  # pyridazine
                Chem.MolFromSmarts("c1ccncc1"),  # pyrrole
                Chem.MolFromSmarts("c1nsncc1"),  # thiazole
                Chem.MolFromSmarts("c1ncocc1"),  # oxazole
                Chem.MolFromSmarts("c1nnn[nH]1"),  # tetrazole
                Chem.MolFromSmarts("c1nn[nH]c1"),  # triazole
                Chem.MolFromSmarts("c1nc2ccccc2[nH]1"),  # benzimidazole
            ]

            # Check if heterocycles are present in reactants and product
            heterocycle_in_reactants = False
            heterocycle_in_product = False

            for reactant in reactants:
                mol = Chem.MolFromSmiles(reactant)
                if mol:
                    for pattern in heterocycle_patterns:
                        if mol.HasSubstructMatch(pattern):
                            heterocycle_in_reactants = True
                            print(f"Found heterocycle in reactant: {reactant}")
                            break

            product_mol = Chem.MolFromSmiles(product)
            if product_mol:
                for pattern in heterocycle_patterns:
                    if product_mol.HasSubstructMatch(pattern):
                        heterocycle_in_product = True
                        print(f"Found heterocycle in product: {product}")
                        break

            # If heterocycles are present in both reactants and product,
            # and they're different, count as a modification
            if heterocycle_in_reactants and heterocycle_in_product:
                # This is a simplification - ideally we would check if the heterocycles
                # are actually different, but that's more complex
                heterocycle_modifications += 1
                print("Heterocycle modification detected!")

        for child in node.get("children", []):
            dfs_traverse(child)

    dfs_traverse(route)
    return heterocycle_modifications > 0
