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
    Detects if the synthesis involves nitration of a pyridine ring.
    """
    has_pyridine_nitration = False

    def dfs_traverse(node):
        nonlocal has_pyridine_nitration

        if node["type"] == "reaction":
            rsmi = node["metadata"]["rsmi"]
            reactants = rsmi.split(">")[0].split(".")
            product = rsmi.split(">")[-1]

            # Check for pyridine in reactants
            pyridine_pattern = Chem.MolFromSmarts("c1ccncc1")
            nitro_pattern = Chem.MolFromSmarts("[N+](=O)[O-]")

            for reactant in reactants:
                reactant_mol = Chem.MolFromSmiles(reactant)
                if reactant_mol and reactant_mol.HasSubstructMatch(pyridine_pattern):
                    # Check if product has pyridine with nitro group
                    product_mol = Chem.MolFromSmiles(product)
                    if (
                        product_mol
                        and product_mol.HasSubstructMatch(pyridine_pattern)
                        and product_mol.HasSubstructMatch(nitro_pattern)
                    ):

                        # Verify reactant doesn't already have nitro group
                        if not reactant_mol.HasSubstructMatch(nitro_pattern):
                            print("Found pyridine nitration")
                            has_pyridine_nitration = True

        for child in node.get("children", []):
            dfs_traverse(child)

    dfs_traverse(route)
    return has_pyridine_nitration
