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
    Detects if the synthetic route involves manipulation of nitrogen-containing heterocycles.
    """
    n_heterocycle_pattern = Chem.MolFromSmarts(
        "c1[n]cc*1"
    )  # Basic N-heterocycle pattern
    has_n_heterocycle_manipulation = False

    def dfs_traverse(node):
        nonlocal has_n_heterocycle_manipulation

        if node["type"] == "reaction":
            if "rsmi" in node.get("metadata", {}):
                rsmi = node["metadata"]["rsmi"]
                reactants = rsmi.split(">")[0].split(".")
                product = rsmi.split(">")[-1]

                try:
                    # Check if both reactants and products contain N-heterocycles
                    reactant_mols = [Chem.MolFromSmiles(r) for r in reactants if r]
                    product_mol = Chem.MolFromSmiles(product)

                    if all(reactant_mols) and product_mol:
                        reactant_has_n_heterocycle = any(
                            m.HasSubstructMatch(n_heterocycle_pattern)
                            for m in reactant_mols
                            if m
                        )
                        product_has_n_heterocycle = product_mol.HasSubstructMatch(
                            n_heterocycle_pattern
                        )

                        if reactant_has_n_heterocycle and product_has_n_heterocycle:
                            # Check if there's a change in the heterocycle
                            main_reactant = max(
                                reactant_mols, key=lambda m: m.GetNumAtoms()
                            )
                            if main_reactant.GetNumAtoms() != product_mol.GetNumAtoms():
                                has_n_heterocycle_manipulation = True
                                print(f"Found N-heterocycle manipulation: {rsmi}")
                except:
                    print(f"Error processing SMILES in N-heterocycle detection: {rsmi}")

        for child in node.get("children", []):
            dfs_traverse(child)

    dfs_traverse(route)
    print(f"N-heterocycle manipulation detected: {has_n_heterocycle_manipulation}")
    return has_n_heterocycle_manipulation
