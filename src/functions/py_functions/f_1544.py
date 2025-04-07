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
    This function detects a sequence involving alcohol protection (as acetate),
    deprotection, and subsequent activation (as mesylate) for nucleophilic substitution.
    """
    acetate_deprotection = False
    alcohol_mesylation = False

    def dfs_traverse(node, depth=0):
        nonlocal acetate_deprotection, alcohol_mesylation

        if node["type"] == "reaction":
            if "rsmi" in node.get("metadata", {}):
                rsmi = node["metadata"]["rsmi"]
                reactants = rsmi.split(">")[0].split(".")
                product = rsmi.split(">")[-1]

                # Check for acetate deprotection
                acetate_pattern = Chem.MolFromSmarts("[C]-[O]-[C](=[O])-[C]")
                alcohol_pattern = Chem.MolFromSmarts("[C]-[O;H1]")

                for reactant in reactants:
                    mol = Chem.MolFromSmiles(reactant)
                    if mol and mol.HasSubstructMatch(acetate_pattern):
                        product_mol = Chem.MolFromSmiles(product)
                        if product_mol and product_mol.HasSubstructMatch(
                            alcohol_pattern
                        ):
                            acetate_deprotection = True
                            print("Acetate deprotection detected at depth", depth)

                # Check for alcohol mesylation
                mesylate_product_pattern = Chem.MolFromSmarts(
                    "[C]-[O]-[S](=[O])(=[O])-[C]"
                )

                for reactant in reactants:
                    mol = Chem.MolFromSmiles(reactant)
                    if mol and mol.HasSubstructMatch(alcohol_pattern):
                        product_mol = Chem.MolFromSmiles(product)
                        if product_mol and product_mol.HasSubstructMatch(
                            mesylate_product_pattern
                        ):
                            alcohol_mesylation = True
                            print("Alcohol mesylation detected at depth", depth)

        # Continue traversing
        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    # Start traversal
    dfs_traverse(route)
    return acetate_deprotection and alcohol_mesylation
