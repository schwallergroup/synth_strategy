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
    Detects if the synthetic route involves reduction of alkyne to alkane linker.
    """
    alkyne_reduction = False

    def dfs_traverse(node):
        nonlocal alkyne_reduction

        if node["type"] == "reaction" and "metadata" in node and "rsmi" in node["metadata"]:
            rsmi = node["metadata"]["rsmi"]
            reactants = rsmi.split(">")[0].split(".")
            product = rsmi.split(">")[-1]

            # Check for alkyne reduction
            for reactant in reactants:
                r_mol = Chem.MolFromSmiles(reactant)
                p_mol = Chem.MolFromSmiles(product)

                if r_mol and p_mol:
                    # Alkyne in reactant
                    if r_mol.HasSubstructMatch(Chem.MolFromSmarts("[c]-[C]#[C]-[c]")):
                        # Corresponding alkane in product
                        if p_mol.HasSubstructMatch(Chem.MolFromSmarts("[c]-[CH2]-[CH2]-[c]")):
                            alkyne_reduction = True
                            print(
                                f"Detected alkyne reduction at depth {node.get('depth', 'unknown')}"
                            )

        # Traverse children
        for child in node.get("children", []):
            dfs_traverse(child)

    dfs_traverse(route)
    return alkyne_reduction
