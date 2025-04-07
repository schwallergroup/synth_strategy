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
    This function detects if the synthesis includes an amide coupling between
    an amine and a carboxylic acid.
    """
    amide_coupling_found = False

    def dfs_traverse(node, depth=0):
        nonlocal amide_coupling_found

        if node["type"] == "reaction":
            if "rsmi" in node.get("metadata", {}):
                rsmi = node["metadata"]["rsmi"]
                reactants = rsmi.split(">")[0].split(".")
                product = rsmi.split(">")[-1]

                # Check if this is an amide coupling reaction
                amine_pattern = Chem.MolFromSmarts("[NH2]")
                acid_pattern = Chem.MolFromSmarts("C(=O)[OH]")
                amide_pattern = Chem.MolFromSmarts("[NH]C(=O)")

                # Check if amine and acid in reactants and amide in product
                product_mol = Chem.MolFromSmiles(product)

                if product_mol and product_mol.HasSubstructMatch(amide_pattern):
                    amine_found = False
                    acid_found = False

                    for reactant in reactants:
                        reactant_mol = Chem.MolFromSmiles(reactant)
                        if reactant_mol:
                            if reactant_mol.HasSubstructMatch(amine_pattern):
                                amine_found = True
                            if reactant_mol.HasSubstructMatch(acid_pattern):
                                acid_found = True

                    if amine_found and acid_found:
                        amide_coupling_found = True
                        print(f"Amide coupling found at depth {depth}")

        # Continue traversing
        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    # Start traversal
    dfs_traverse(route)
    return amide_coupling_found
