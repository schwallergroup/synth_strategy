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
    This function detects if the synthetic route involves multiple functional group
    transformations on aromatic rings.
    """
    aromatic_transformations = 0

    def dfs_traverse(node):
        nonlocal aromatic_transformations

        if node["type"] == "reaction" and "metadata" in node and "rsmi" in node["metadata"]:
            rsmi = node["metadata"]["rsmi"]
            reactants = rsmi.split(">")[0].split(".")
            product = rsmi.split(">")[-1]

            reactant_mols = [Chem.MolFromSmiles(r) for r in reactants if r]
            product_mol = Chem.MolFromSmiles(product) if product else None

            if product_mol and reactant_mols:
                # Check for aromatic rings in reactants
                aromatic_pattern = Chem.MolFromSmarts("c1ccccc1")
                has_aromatic_reactants = any(
                    mol.HasSubstructMatch(aromatic_pattern) for mol in reactant_mols if mol
                )

                # Check for aromatic rings in product
                has_aromatic_product = (
                    product_mol.HasSubstructMatch(aromatic_pattern) if product_mol else False
                )

                if has_aromatic_reactants and has_aromatic_product:
                    # Check for functional group transformations on aromatic rings

                    # Nitro reduction
                    nitro_pattern = Chem.MolFromSmarts("c-[NX3+](=[OX1])[OX1-]")
                    amine_pattern = Chem.MolFromSmarts("c-[NX3;H2]")

                    has_nitro = any(
                        mol.HasSubstructMatch(nitro_pattern) for mol in reactant_mols if mol
                    )
                    has_amine = (
                        product_mol.HasSubstructMatch(amine_pattern) if product_mol else False
                    )

                    if has_nitro and has_amine:
                        aromatic_transformations += 1
                        print("Aromatic nitro reduction detected")

                    # Amide formation on aromatic ring
                    aromatic_acid_pattern = Chem.MolFromSmarts("c-[CX3](=[OX1])[OX2H]")
                    aromatic_amide_pattern = Chem.MolFromSmarts("c-[CX3](=[OX1])[NX3]")

                    has_aromatic_acid = any(
                        mol.HasSubstructMatch(aromatic_acid_pattern) for mol in reactant_mols if mol
                    )
                    has_aromatic_amide = (
                        product_mol.HasSubstructMatch(aromatic_amide_pattern)
                        if product_mol
                        else False
                    )

                    if has_aromatic_acid and has_aromatic_amide:
                        aromatic_transformations += 1
                        print("Aromatic amide formation detected")

        for child in node.get("children", []):
            dfs_traverse(child)

    dfs_traverse(route)
    return aromatic_transformations >= 2
