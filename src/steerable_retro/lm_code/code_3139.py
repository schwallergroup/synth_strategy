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
    This function detects if the synthetic route involves a linker extension strategy
    where a multi-carbon linker is introduced between an aromatic ring and an amine.
    """
    linker_extension = False

    def dfs_traverse(node):
        nonlocal linker_extension

        if node["type"] == "reaction" and "metadata" in node and "rsmi" in node["metadata"]:
            rsmi = node["metadata"]["rsmi"]
            reactants = rsmi.split(">")[0].split(".")
            product = rsmi.split(">")[-1]

            try:
                prod_mol = Chem.MolFromSmiles(product)
                # Look for pattern where aromatic ring is connected to amine via 3+ atom linker
                linker_pattern = Chem.MolFromSmarts("c[O][CH2][CH2][CH2][#7]")

                if prod_mol and prod_mol.HasSubstructMatch(linker_pattern):
                    # Check if this pattern wasn't present in reactants
                    pattern_in_reactants = False
                    for r in reactants:
                        r_mol = Chem.MolFromSmiles(r)
                        if r_mol and r_mol.HasSubstructMatch(linker_pattern):
                            pattern_in_reactants = True
                            break

                    if not pattern_in_reactants:
                        linker_extension = True
                        print(
                            "Detected linker extension strategy connecting aromatic ring to amine"
                        )
            except:
                pass

        for child in node.get("children", []):
            dfs_traverse(child)

    dfs_traverse(route)
    return linker_extension
