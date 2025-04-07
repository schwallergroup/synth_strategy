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
    Detects if the synthesis route involves aromatic fluorination.
    """
    found_fluorination = False

    def dfs_traverse(node):
        nonlocal found_fluorination

        if node["type"] == "reaction" and "metadata" in node and "rsmi" in node["metadata"]:
            rsmi = node["metadata"]["rsmi"]
            reactants = rsmi.split(">")[0].split(".")
            product = rsmi.split(">")[-1]

            # Check for fluorination pattern
            aryl_fluoride_pattern = Chem.MolFromSmarts("[c]-[F]")

            reactant_mols = [Chem.MolFromSmiles(r) for r in reactants if r]
            product_mol = Chem.MolFromSmiles(product) if product else None

            # Count fluorine atoms in reactants and product
            if product_mol:
                product_f_count = len(product_mol.GetSubstructMatches(aryl_fluoride_pattern))
                reactant_f_count = sum(
                    len(mol.GetSubstructMatches(aryl_fluoride_pattern))
                    for mol in reactant_mols
                    if mol
                )

                if product_f_count > reactant_f_count:
                    found_fluorination = True
                    print("Found aromatic fluorination step")

        for child in node.get("children", []):
            dfs_traverse(child)

    dfs_traverse(route)
    return found_fluorination
