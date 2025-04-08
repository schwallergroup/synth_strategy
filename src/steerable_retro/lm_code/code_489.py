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
    Detects if the synthesis route involves an amide bond formation from
    a carboxylic acid and an amine.
    """
    amide_formation_found = False

    def dfs_traverse(node):
        nonlocal amide_formation_found

        if node.get("type") == "reaction" and "metadata" in node and "rsmi" in node["metadata"]:
            rsmi = node["metadata"]["rsmi"]
            reactants_smiles = rsmi.split(">")[0]
            product_smiles = rsmi.split(">")[-1]

            # Check for carboxylic acid and amine in reactants, amide in product
            carboxylic_acid_pattern = Chem.MolFromSmarts("[C](=O)[OH]")
            amine_pattern = Chem.MolFromSmarts("[N;!$(NC=O)]")
            amide_pattern = Chem.MolFromSmarts("[C](=O)[N]")

            if carboxylic_acid_pattern and amine_pattern and amide_pattern:
                reactants_mols = [Chem.MolFromSmiles(r) for r in reactants_smiles.split(".") if r]
                product_mol = Chem.MolFromSmiles(product_smiles)

                if product_mol and any(reactants_mols):
                    has_acid = any(
                        r.HasSubstructMatch(carboxylic_acid_pattern) for r in reactants_mols if r
                    )
                    has_amine = any(r.HasSubstructMatch(amine_pattern) for r in reactants_mols if r)
                    has_amide = product_mol.HasSubstructMatch(amide_pattern)

                    if has_acid and has_amine and has_amide:
                        amide_formation_found = True
                        print("Amide bond formation detected")

        for child in node.get("children", []):
            dfs_traverse(child)

    dfs_traverse(route)
    return amide_formation_found
