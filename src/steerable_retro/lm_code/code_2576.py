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
    Detects if the synthesis route involves SNAr reaction with aryl chloride and aromatic amine.
    """
    found_snar = False

    def dfs_traverse(node):
        nonlocal found_snar

        if node["type"] == "reaction" and "metadata" in node and "rsmi" in node["metadata"]:
            rsmi = node["metadata"]["rsmi"]
            reactants = rsmi.split(">")[0].split(".")
            product = rsmi.split(">")[-1]

            # Check for SNAr pattern
            aryl_chloride_pattern = Chem.MolFromSmarts("[c]-[Cl]")
            aromatic_amine_pattern = Chem.MolFromSmarts("[NH2]-[c]")
            diarylamine_pattern = Chem.MolFromSmarts("[c]-[NH]-[c]")

            reactant_mols = [Chem.MolFromSmiles(r) for r in reactants if r]
            product_mol = Chem.MolFromSmiles(product) if product else None

            if (
                product_mol
                and product_mol.HasSubstructMatch(diarylamine_pattern)
                and any(
                    mol and mol.HasSubstructMatch(aryl_chloride_pattern) for mol in reactant_mols
                )
                and any(
                    mol and mol.HasSubstructMatch(aromatic_amine_pattern) for mol in reactant_mols
                )
            ):
                found_snar = True
                print("Found SNAr reaction step")

        for child in node.get("children", []):
            dfs_traverse(child)

    dfs_traverse(route)
    return found_snar
