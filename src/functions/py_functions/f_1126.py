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
    This function detects if the synthesis includes hydrazine formation.
    """
    hydrazine_formation_found = False

    def dfs_traverse(node):
        nonlocal hydrazine_formation_found

        if node["type"] == "reaction":
            if "metadata" in node and "rsmi" in node["metadata"]:
                rsmi = node["metadata"]["rsmi"]
                reactants = rsmi.split(">")[0].split(".")
                product = rsmi.split(">")[-1]

                # Check for hydrazine in reactants
                hydrazine_in_reactants = False
                for reactant in reactants:
                    reactant_mol = Chem.MolFromSmiles(reactant)
                    if reactant_mol:
                        hydrazine_pattern = Chem.MolFromSmarts("[NH2][NH2]")
                        if reactant_mol.HasSubstructMatch(hydrazine_pattern):
                            hydrazine_in_reactants = True
                            break

                # Check for N-N bond in product
                if hydrazine_in_reactants:
                    product_mol = Chem.MolFromSmiles(product)
                    if product_mol:
                        nn_bond_pattern = Chem.MolFromSmarts("[N]-[N]")
                        if product_mol.HasSubstructMatch(nn_bond_pattern):
                            hydrazine_formation_found = True
                            print("Detected hydrazine formation")

        for child in node.get("children", []):
            dfs_traverse(child)

    dfs_traverse(route)
    return hydrazine_formation_found
