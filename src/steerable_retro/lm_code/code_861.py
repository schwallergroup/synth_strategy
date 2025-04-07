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
    Detects a synthetic route that includes nitration of an aromatic system.
    """
    has_aromatic_nitration = False

    def dfs_traverse(node):
        nonlocal has_aromatic_nitration

        if node["type"] == "reaction":
            rsmi = node["metadata"]["rsmi"]
            reactants_smiles = rsmi.split(">")[0].split(".")
            product_smiles = rsmi.split(">")[-1]

            reactant_mols = [Chem.MolFromSmiles(r) for r in reactants_smiles]
            product_mol = Chem.MolFromSmiles(product_smiles)

            # Check for nitric acid or nitronium ion
            nitric_acid_pattern = Chem.MolFromSmarts("O[N+](=O)[O-]")

            # Check for aromatic system in reactants and nitro-aromatic in product
            aromatic_pattern = Chem.MolFromSmarts("c")
            nitro_aromatic_pattern = Chem.MolFromSmarts("c[N+](=O)[O-]")

            has_nitric_acid = any(
                r and r.HasSubstructMatch(nitric_acid_pattern) for r in reactant_mols
            )
            has_aromatic = any(r and r.HasSubstructMatch(aromatic_pattern) for r in reactant_mols)
            has_nitro_aromatic = product_mol and product_mol.HasSubstructMatch(
                nitro_aromatic_pattern
            )

            if (
                (has_nitric_acid or "NO2" in "".join(reactants_smiles))
                and has_aromatic
                and has_nitro_aromatic
            ):
                has_aromatic_nitration = True
                print("Found aromatic nitration")

        for child in node.get("children", []):
            dfs_traverse(child)

    dfs_traverse(route)
    print(f"Aromatic nitration strategy: {has_aromatic_nitration}")
    return has_aromatic_nitration
