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
    This function detects phthalimide deprotection to reveal a primary amine.
    """
    phthalimide_deprotection_detected = False

    def dfs_traverse(node):
        nonlocal phthalimide_deprotection_detected

        if node["type"] == "reaction":
            if "rsmi" in node.get("metadata", {}):
                rsmi = node["metadata"]["rsmi"]
                reactants = rsmi.split(">")[0].split(".")
                product = rsmi.split(">")[-1]

                # Check for phthalimide pattern in reactants
                phthalimide_pattern = Chem.MolFromSmarts("N1C(=O)c2ccccc2C1=O")
                amine_pattern = Chem.MolFromSmarts("[NH2]")

                reactant_mols = [Chem.MolFromSmiles(r) for r in reactants if r]
                product_mol = Chem.MolFromSmiles(product) if product else None

                if product_mol and any(
                    m and m.HasSubstructMatch(phthalimide_pattern)
                    for m in reactant_mols
                ):
                    if product_mol.HasSubstructMatch(amine_pattern):
                        print("Phthalimide deprotection detected")
                        phthalimide_deprotection_detected = True

        for child in node.get("children", []):
            dfs_traverse(child)

    dfs_traverse(route)
    return phthalimide_deprotection_detected
