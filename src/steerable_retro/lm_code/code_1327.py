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
    Detects if the synthesis involves a nitro reduction to amine.
    """
    has_nitro_reduction = False

    def dfs_traverse(node):
        nonlocal has_nitro_reduction

        if node["type"] == "reaction":
            rsmi = node["metadata"]["rsmi"]
            reactants = rsmi.split(">")[0].split(".")
            product = rsmi.split(">")[-1]

            # Check for nitro in reactant and amine in product
            for reactant in reactants:
                reactant_mol = Chem.MolFromSmiles(reactant)
                product_mol = Chem.MolFromSmiles(product)

                if reactant_mol and product_mol:
                    nitro_pattern = Chem.MolFromSmarts("[#6][N+](=[O])[O-]")
                    amine_pattern = Chem.MolFromSmarts("[#6][NH2]")

                    has_nitro = reactant_mol.HasSubstructMatch(nitro_pattern)
                    has_amine = product_mol.HasSubstructMatch(amine_pattern)

                    if has_nitro and has_amine:
                        print(
                            f"Found nitro reduction at depth {node.get('metadata', {}).get('depth', 'unknown')}"
                        )
                        has_nitro_reduction = True

        for child in node.get("children", []):
            dfs_traverse(child)

    dfs_traverse(route)
    return has_nitro_reduction
