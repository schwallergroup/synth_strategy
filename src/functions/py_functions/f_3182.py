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
    Detects if the synthesis route involves amine protection/deprotection strategies.
    """
    has_protection = False

    def dfs_traverse(node):
        nonlocal has_protection

        if node["type"] == "reaction":
            # Extract reactants and product
            rsmi = node["metadata"]["rsmi"]
            reactants_smiles = rsmi.split(">")[0].split(".")
            product_smiles = rsmi.split(">")[-1]

            # Check for Boc protection
            product_mol = Chem.MolFromSmiles(product_smiles)
            if product_mol is not None:
                # Boc protected amine pattern
                boc_pattern = Chem.MolFromSmarts("[#7]C(=O)OC(C)(C)C")
                if product_mol.HasSubstructMatch(boc_pattern):
                    # Check if reactants had primary/secondary amine
                    reactant_mols = [Chem.MolFromSmiles(r) for r in reactants_smiles]
                    amine_pattern = Chem.MolFromSmarts("[#7;H1,H2]")
                    if any(
                        mol is not None and mol.HasSubstructMatch(amine_pattern)
                        for mol in reactant_mols
                    ):
                        has_protection = True
                        print("Detected Boc protection of amine")

        # Traverse children
        for child in node.get("children", []):
            dfs_traverse(child)

    # Start traversal from root
    dfs_traverse(route)
    return has_protection
