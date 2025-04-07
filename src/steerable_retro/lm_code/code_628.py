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
    This function detects a synthetic strategy centered around pyrazine derivatives,
    where pyrazine is present throughout the synthesis and undergoes multiple
    functionalization steps.
    """
    pyrazine_count = 0
    has_pyrazine_functionalization = False

    def dfs_traverse(node):
        nonlocal pyrazine_count, has_pyrazine_functionalization

        if node["type"] == "mol" and "smiles" in node:
            try:
                mol = Chem.MolFromSmiles(node["smiles"])
                if mol:
                    pyrazine_pattern = Chem.MolFromSmarts("c1ncccn1")
                    if mol.HasSubstructMatch(pyrazine_pattern):
                        pyrazine_count += 1
            except:
                pass

        elif node["type"] == "reaction" and "metadata" in node and "rsmi" in node["metadata"]:
            rsmi = node["metadata"]["rsmi"]
            reactants_smiles = rsmi.split(">")[0].split(".")
            product_smiles = rsmi.split(">")[-1]

            # Check for pyrazine functionalization
            pyrazine_pattern = Chem.MolFromSmarts("c1ncccn1")

            try:
                product_mol = Chem.MolFromSmiles(product_smiles)
                if product_mol and product_mol.HasSubstructMatch(pyrazine_pattern):
                    for reactant in reactants_smiles:
                        reactant_mol = Chem.MolFromSmiles(reactant)
                        if reactant_mol and reactant_mol.HasSubstructMatch(pyrazine_pattern):
                            # If both reactant and product have pyrazine, it's a functionalization
                            has_pyrazine_functionalization = True
                            print("Detected pyrazine functionalization")
                            break
            except:
                pass

        # Continue traversing
        for child in node.get("children", []):
            dfs_traverse(child)

    # Start traversal from the root
    dfs_traverse(route)

    # Return True if pyrazine is present in multiple steps and undergoes functionalization
    return pyrazine_count >= 3 and has_pyrazine_functionalization
