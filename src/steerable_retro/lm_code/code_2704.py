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
    This function detects if the synthesis route includes aromatization of a heterocycle.
    """
    aromatization_found = False

    def dfs_traverse(node):
        nonlocal aromatization_found

        if node["type"] == "reaction" and "metadata" in node and "rsmi" in node["metadata"]:
            rsmi = node["metadata"]["rsmi"]
            reactants = rsmi.split(">")[0].split(".")
            product = rsmi.split(">")[-1]

            # Check for non-aromatic heterocycle in reactants and aromatic heterocycle in product
            non_aromatic_n_ring = Chem.MolFromSmarts("[#6]1~[#6]~[#6]~[#6]~[#7]~[#6]1")
            aromatic_n_ring = Chem.MolFromSmarts("[#6]1:[#6]:[#6]:[#6]:[#7]:[#6]:1")

            try:
                for reactant in reactants:
                    reactant_mol = Chem.MolFromSmiles(reactant)
                    if reactant_mol and reactant_mol.HasSubstructMatch(non_aromatic_n_ring):
                        product_mol = Chem.MolFromSmiles(product)
                        if product_mol and product_mol.HasSubstructMatch(aromatic_n_ring):
                            print("Heterocycle aromatization detected")
                            aromatization_found = True
                            break
            except:
                pass

        # Traverse children
        for child in node.get("children", []):
            dfs_traverse(child)

    # Start traversal from the root
    dfs_traverse(route)
    return aromatization_found
