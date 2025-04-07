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
    Detects if the synthesis involves O-methylation of a hydroxyl group.
    """
    result = False

    def dfs_traverse(node):
        nonlocal result

        if node["type"] == "reaction":
            # Extract reactants and product
            rsmi = node["metadata"]["rsmi"]
            reactants_smiles = rsmi.split(">")[0].split(".")
            product_smiles = rsmi.split(">")[-1]

            # Check for hydroxyl in reactants
            hydroxyl_pattern = Chem.MolFromSmarts("[OH]")

            # Check for methoxy in product
            methoxy_pattern = Chem.MolFromSmarts("[O][CH3]")

            has_hydroxyl = False
            for reactant in reactants_smiles:
                try:
                    mol = Chem.MolFromSmiles(reactant)
                    if mol and mol.HasSubstructMatch(hydroxyl_pattern):
                        has_hydroxyl = True
                        break
                except:
                    continue

            if has_hydroxyl:
                try:
                    product_mol = Chem.MolFromSmiles(product_smiles)
                    if product_mol and product_mol.HasSubstructMatch(methoxy_pattern):
                        # Check if methoxy count in product is greater than in reactants
                        methoxy_in_reactants = 0
                        for reactant in reactants_smiles:
                            try:
                                r_mol = Chem.MolFromSmiles(reactant)
                                if r_mol:
                                    methoxy_in_reactants += len(
                                        r_mol.GetSubstructMatches(methoxy_pattern)
                                    )
                            except:
                                continue

                        methoxy_in_product = len(
                            product_mol.GetSubstructMatches(methoxy_pattern)
                        )

                        if methoxy_in_product > methoxy_in_reactants:
                            print("Detected O-methylation of hydroxyl group")
                            result = True
                except:
                    pass

        # Continue traversing
        for child in node.get("children", []):
            dfs_traverse(child)

    dfs_traverse(route)
    return result
