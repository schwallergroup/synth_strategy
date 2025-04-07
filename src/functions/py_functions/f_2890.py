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
    This function detects the formation of a diaryl ether linkage.
    It looks for a reaction where a C-O-C bond is formed between two aromatic rings.
    """
    diaryl_ether_formed = False

    def dfs_traverse(node):
        nonlocal diaryl_ether_formed

        if node["type"] == "reaction" and "rsmi" in node.get("metadata", {}):
            rsmi = node["metadata"]["rsmi"]
            reactants = rsmi.split(">")[0].split(".")
            product = rsmi.split(">")[-1]

            # Check for diaryl ether formation
            diaryl_ether_pattern = Chem.MolFromSmarts("c-[O]-c")
            phenol_pattern = Chem.MolFromSmarts("c-[OH]")
            fluoro_arene_pattern = Chem.MolFromSmarts("c-[F]")

            reactant_mols = [Chem.MolFromSmiles(r) for r in reactants if r]
            product_mol = Chem.MolFromSmiles(product) if product else None

            if product_mol and product_mol.HasSubstructMatch(diaryl_ether_pattern):
                # Check if reactants contain phenol and fluoro-arene
                has_phenol = any(
                    m and m.HasSubstructMatch(phenol_pattern) for m in reactant_mols
                )
                has_fluoro_arene = any(
                    m and m.HasSubstructMatch(fluoro_arene_pattern)
                    for m in reactant_mols
                )

                if has_phenol and has_fluoro_arene:
                    diaryl_ether_formed = True
                    print("Found diaryl ether formation from phenol and fluoro-arene")

        for child in node.get("children", []):
            dfs_traverse(child)

    dfs_traverse(route)
    return diaryl_ether_formed
