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
    Detects formation of heterocyclic systems like benzoxazinone.
    """
    heterocycle_formation = False

    def dfs_traverse(node, depth=0):
        nonlocal heterocycle_formation

        if node["type"] == "reaction":
            if "rsmi" in node.get("metadata", {}):
                rsmi = node["metadata"]["rsmi"]
                reactants_smiles = rsmi.split(">")[0].split(".")
                product_smiles = rsmi.split(">")[-1]

                try:
                    reactant_mols = [Chem.MolFromSmiles(r) for r in reactants_smiles]
                    product_mol = Chem.MolFromSmiles(product_smiles)

                    if (
                        all(mol is not None for mol in reactant_mols)
                        and product_mol is not None
                    ):
                        # Check for benzoxazinone formation
                        benzoxazinone_pattern = Chem.MolFromSmarts(
                            "c1ccc2c(c1)OCC(=O)N2"
                        )

                        has_benzoxazinone = product_mol.HasSubstructMatch(
                            benzoxazinone_pattern
                        )
                        reactants_have_benzoxazinone = any(
                            mol.HasSubstructMatch(benzoxazinone_pattern)
                            for mol in reactant_mols
                        )

                        if has_benzoxazinone and not reactants_have_benzoxazinone:
                            heterocycle_formation = True
                            print(f"Heterocycle formation detected at depth {depth}")
                except Exception as e:
                    print(f"Error processing reaction: {e}")

        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    dfs_traverse(route)

    return heterocycle_formation
