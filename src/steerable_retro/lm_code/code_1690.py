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
    This function detects if a fluorinated aromatic ring is introduced via a coupling reaction.
    """
    fluorinated_coupling = False

    def dfs_traverse(node):
        nonlocal fluorinated_coupling

        if node["type"] == "reaction" and "metadata" in node and "rsmi" in node["metadata"]:
            rsmi = node["metadata"]["rsmi"]
            reactants = rsmi.split(">")[0].split(".")
            product = rsmi.split(">")[-1]

            # Check for fluorinated aromatic in reactants
            fluoro_aromatic_pattern = Chem.MolFromSmarts("[c][F]")

            # Check for coupling reaction (biaryl formation)
            biaryl_pattern = Chem.MolFromSmarts("[c]!@[c]")

            try:
                has_fluoro_aromatic = False
                for reactant in reactants:
                    mol = Chem.MolFromSmiles(reactant)
                    if mol and mol.HasSubstructMatch(fluoro_aromatic_pattern):
                        has_fluoro_aromatic = True
                        break

                product_mol = Chem.MolFromSmiles(product)
                has_biaryl = product_mol and product_mol.HasSubstructMatch(biaryl_pattern)
                has_fluoro_product = product_mol and product_mol.HasSubstructMatch(
                    fluoro_aromatic_pattern
                )

                if has_fluoro_aromatic and has_biaryl and has_fluoro_product:
                    print("Fluorinated aromatic introduced via coupling reaction")
                    fluorinated_coupling = True
            except:
                pass

        for child in node.get("children", []):
            dfs_traverse(child)

    dfs_traverse(route)
    return fluorinated_coupling
