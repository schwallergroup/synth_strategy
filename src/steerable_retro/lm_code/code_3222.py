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
    This function detects if the final product contains both fluorine and iodine on aromatic rings.
    """
    has_fluoroiodo_aromatic = False

    def dfs_traverse(node):
        nonlocal has_fluoroiodo_aromatic

        if node["type"] == "mol" and node.get("in_stock", False) == False:
            # This is likely the final product
            try:
                mol = Chem.MolFromSmiles(node["smiles"])
                if mol:
                    # Check for both F and I on aromatic rings
                    fluoro_aromatic = mol.HasSubstructMatch(Chem.MolFromSmarts("[c]-[F]"))
                    iodo_aromatic = mol.HasSubstructMatch(Chem.MolFromSmarts("[c]-[I]"))

                    if fluoro_aromatic and iodo_aromatic:
                        print(
                            "Detected both fluorine and iodine on aromatic rings in molecule:",
                            node["smiles"],
                        )
                        has_fluoroiodo_aromatic = True
            except:
                pass

        for child in node.get("children", []):
            dfs_traverse(child)

    dfs_traverse(route)
    return has_fluoroiodo_aromatic
