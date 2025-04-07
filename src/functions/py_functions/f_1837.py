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
    Detects if the synthesis route includes aromatic bromination.
    Looks for a reaction where a bromine atom is added to an aromatic carbon.
    """
    bromination_found = False

    def dfs_traverse(node):
        nonlocal bromination_found

        if (
            node["type"] == "reaction"
            and "metadata" in node
            and "rsmi" in node["metadata"]
        ):
            rsmi = node["metadata"]["rsmi"]
            reactants = rsmi.split(">")[0].split(".")
            product = rsmi.split(">")[-1]

            # Check if product has a bromine on an aromatic carbon that wasn't in reactants
            product_mol = Chem.MolFromSmiles(product)
            if product_mol:
                # Look for aromatic carbon with bromine
                patt = Chem.MolFromSmarts("[c][Br]")
                if product_mol.HasSubstructMatch(patt):
                    # Check if this is new (not in reactants)
                    br_in_reactants = False
                    for reactant in reactants:
                        reactant_mol = Chem.MolFromSmiles(reactant)
                        if reactant_mol and reactant_mol.HasSubstructMatch(patt):
                            br_in_reactants = True
                            break

                    if not br_in_reactants:
                        print(f"Aromatic bromination detected: {rsmi}")
                        bromination_found = True

        for child in node.get("children", []):
            dfs_traverse(child)

    dfs_traverse(route)
    return bromination_found
