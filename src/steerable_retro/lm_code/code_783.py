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
    This function detects SNAr reactions with fluorine as a leaving group.
    """
    snar_with_f_found = False

    def dfs_traverse(node):
        nonlocal snar_with_f_found

        if node["type"] == "reaction":
            # Extract reactants and product from reaction SMILES
            rsmi = node["metadata"]["rsmi"]
            reactants_smiles = rsmi.split(">")[0].split(".")
            product_smiles = rsmi.split(">")[-1]

            # Check for aromatic-F bond in reactants
            f_pattern = Chem.MolFromSmarts("c-[F]")
            f_reactant = None

            for reactant in reactants_smiles:
                try:
                    mol = Chem.MolFromSmiles(reactant)
                    if mol and mol.HasSubstructMatch(f_pattern):
                        f_reactant = reactant
                        break
                except:
                    continue

            if f_reactant:
                # Check if the F is replaced in the product (typically by N or O)
                # This is a simplification - a more robust implementation would track atom mappings
                product_mol = Chem.MolFromSmiles(product_smiles)
                if product_mol:
                    # Look for new C-N or C-O bond where F was attached
                    # This is a heuristic approach
                    if (
                        len(reactants_smiles) > 1
                        and ("N" in product_smiles or "O" in product_smiles)
                        and not product_mol.HasSubstructMatch(f_pattern)
                    ):
                        print("SNAr with fluorine leaving group detected")
                        snar_with_f_found = True

        # Traverse children
        for child in node.get("children", []):
            dfs_traverse(child)

    # Start traversal from root
    dfs_traverse(route)
    return snar_with_f_found
