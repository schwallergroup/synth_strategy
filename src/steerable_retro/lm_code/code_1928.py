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
    This function detects benzylic functionalization, specifically bromination
    of a benzylic position.
    """
    has_benzylic_bromination = False

    def dfs_traverse(node):
        nonlocal has_benzylic_bromination

        if node["type"] == "reaction":
            if "rsmi" in node.get("metadata", {}):
                rsmi = node["metadata"]["rsmi"]
                reactants = rsmi.split(">")[0].split(".")
                product = rsmi.split(">")[-1]

                # Check for benzylic methyl in reactants
                benzylic_methyl_pattern = Chem.MolFromSmarts("c[CH3]")
                benzylic_bromide_pattern = Chem.MolFromSmarts("c[CH2][Br]")

                has_benzylic_methyl = False
                for reactant in reactants:
                    mol = Chem.MolFromSmiles(reactant)
                    if mol and mol.HasSubstructMatch(benzylic_methyl_pattern):
                        has_benzylic_methyl = True
                        break

                # Check if product has benzylic bromide
                product_mol = Chem.MolFromSmiles(product)
                if (
                    has_benzylic_methyl
                    and product_mol
                    and product_mol.HasSubstructMatch(benzylic_bromide_pattern)
                ):
                    has_benzylic_bromination = True
                    print("Detected benzylic bromination")

        # Traverse children
        for child in node.get("children", []):
            dfs_traverse(child)

    # Start traversal from root
    dfs_traverse(route)
    return has_benzylic_bromination
