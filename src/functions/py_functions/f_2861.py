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
    This function detects an early-stage heterocycle transformation from pyran to pyridine.
    """
    found_transformation = False

    def dfs_traverse(node, depth=0):
        nonlocal found_transformation

        if node["type"] == "reaction" and depth >= 3:  # Early stage (high depth)
            if "rsmi" in node.get("metadata", {}):
                rsmi = node["metadata"]["rsmi"]
                reactants = rsmi.split(">")[0].split(".")
                product = rsmi.split(">")[-1]

                # Check for pyran pattern in reactants
                pyran_pattern = Chem.MolFromSmarts("o1ccccc1=O")
                # Check for pyridine pattern in product
                pyridine_pattern = Chem.MolFromSmarts("n1ccccc1=O")

                for reactant in reactants:
                    reactant_mol = Chem.MolFromSmiles(reactant)
                    if reactant_mol and reactant_mol.HasSubstructMatch(pyran_pattern):
                        product_mol = Chem.MolFromSmiles(product)
                        if product_mol and product_mol.HasSubstructMatch(
                            pyridine_pattern
                        ):
                            print(
                                f"Found pyran to pyridine transformation at depth {depth}"
                            )
                            found_transformation = True

        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    dfs_traverse(route)
    return found_transformation
