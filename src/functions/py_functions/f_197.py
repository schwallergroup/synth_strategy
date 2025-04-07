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
    Detects if the synthesis involves formation of a heterocyclic ring system.
    """
    ring_formation_found = False

    def dfs_traverse(node):
        nonlocal ring_formation_found

        if (
            node["type"] == "reaction"
            and "metadata" in node
            and "rsmi" in node["metadata"]
        ):
            rsmi = node["metadata"]["rsmi"]
            reactants = rsmi.split(">")[0].split(".")
            product = rsmi.split(">")[-1]

            # Check for ring formation by comparing ring count
            reactant_mol = Chem.MolFromSmiles(reactants[0])
            product_mol = Chem.MolFromSmiles(product)

            if reactant_mol and product_mol:
                reactant_rings = Chem.GetSSSR(reactant_mol)
                product_rings = Chem.GetSSSR(product_mol)

                if len(product_rings) > len(reactant_rings):
                    # Check if new ring contains a nitrogen (heterocycle)
                    nitrogen_pattern = Chem.MolFromSmarts("[n]")
                    if product_mol.HasSubstructMatch(nitrogen_pattern):
                        ring_formation_found = True
                        print(
                            f"Found heterocyclic ring formation at depth {node.get('metadata', {}).get('ID', '')}"
                        )

        for child in node.get("children", []):
            dfs_traverse(child)

    dfs_traverse(route)
    return ring_formation_found
