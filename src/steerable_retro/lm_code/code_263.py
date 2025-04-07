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
    This function detects amide formation via coupling of carboxylic acid with amine.
    """
    amide_coupling_detected = False

    def dfs_traverse(node):
        nonlocal amide_coupling_detected

        if node["type"] == "reaction":
            # Extract reactants and product
            rsmi = node["metadata"]["rsmi"]
            reactants_smiles = rsmi.split(">")[0].split(".")
            product_smiles = rsmi.split(">")[-1]

            # Check for carboxylic acid and amine in reactants, amide in product
            acid_pattern = Chem.MolFromSmarts("[C;$(C=O)][OH]")
            amine_pattern = Chem.MolFromSmarts("[N;!$(N=*);!$(NC=O)]")
            amide_pattern = Chem.MolFromSmarts("[C;$(C=O)][N;!$(N=*)]")

            reactants_mols = [Chem.MolFromSmiles(r) for r in reactants_smiles]
            product_mol = Chem.MolFromSmiles(product_smiles)

            if product_mol and amide_pattern:
                has_acid = any(r and r.HasSubstructMatch(acid_pattern) for r in reactants_mols if r)
                has_amine = any(
                    r and r.HasSubstructMatch(amine_pattern) for r in reactants_mols if r
                )
                product_has_amide = product_mol.HasSubstructMatch(amide_pattern)

                if has_acid and has_amine and product_has_amide:
                    amide_coupling_detected = True
                    print("Amide coupling with amine detected")

        # Traverse children
        for child in node.get("children", []):
            dfs_traverse(child)

    # Start traversal
    dfs_traverse(route)

    return amide_coupling_detected
