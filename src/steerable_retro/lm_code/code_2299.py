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
    Detects if the route contains a convergent ether formation step where
    a phenol and an alcohol are combined to form an ether.
    """
    found_ether_formation = False

    def dfs_traverse(node, depth=0):
        nonlocal found_ether_formation

        if node["type"] == "reaction":
            if "rsmi" in node["metadata"]:
                rsmi = node["metadata"]["rsmi"]
                reactants = rsmi.split(">")[0].split(".")
                product = rsmi.split(">")[-1]

                # Check for phenol and alcohol in reactants
                phenol_pattern = Chem.MolFromSmarts("c[OH]")
                alcohol_pattern = Chem.MolFromSmarts("C[OH]")
                ether_pattern = Chem.MolFromSmarts("c[O]C")

                has_phenol = False
                has_alcohol = False

                for r in reactants:
                    try:
                        mol = Chem.MolFromSmiles(r)
                        if mol:
                            if mol.HasSubstructMatch(phenol_pattern):
                                has_phenol = True
                            if mol.HasSubstructMatch(alcohol_pattern):
                                has_alcohol = True
                    except:
                        continue

                try:
                    product_mol = Chem.MolFromSmiles(product)
                    if (
                        has_phenol
                        and has_alcohol
                        and product_mol
                        and product_mol.HasSubstructMatch(ether_pattern)
                    ):
                        found_ether_formation = True
                        print(f"Found convergent ether formation at depth {depth}")
                except:
                    pass

        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    dfs_traverse(route)
    return found_ether_formation
