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
    Detects a synthesis strategy that uses alkyl bromide linkers (particularly 1,3-dibromopropane)
    to connect fragments through ether formations.
    """
    # Track key features
    alkyl_bromide_linker_used = False
    ether_formations = 0

    # Define SMARTS patterns
    alkyl_bromide_pattern = Chem.MolFromSmarts("[Br][CH2][CH2][CH2][Br]")
    mono_bromide_pattern = Chem.MolFromSmarts("[Br][CH2][CH2][CH2][O]")
    ether_pattern = Chem.MolFromSmarts("[#6]-[O]-[#6]")

    def dfs_traverse(node, depth=0):
        nonlocal alkyl_bromide_linker_used, ether_formations

        if node["type"] == "reaction":
            # Extract reactants and product
            rsmi = node["metadata"]["rsmi"]
            reactants_smiles = rsmi.split(">")[0].split(".")
            product_smiles = rsmi.split(">")[-1]

            # Check for alkyl bromide linker
            reactant_mols = [Chem.MolFromSmiles(r) for r in reactants_smiles]
            product_mol = Chem.MolFromSmiles(product_smiles)

            for reactant in reactant_mols:
                if reactant and reactant.HasSubstructMatch(alkyl_bromide_pattern):
                    alkyl_bromide_linker_used = True
                    print(f"1,3-dibromopropane linker detected at depth {depth}")

            # Check for mono-bromide intermediate (partially coupled linker)
            for reactant in reactant_mols:
                if reactant and reactant.HasSubstructMatch(mono_bromide_pattern):
                    print(f"Partially coupled linker detected at depth {depth}")

            # Count ethers in reactants and product
            reactant_ethers = sum(
                len(mol.GetSubstructMatches(ether_pattern)) for mol in reactant_mols if mol
            )
            product_ethers = (
                len(product_mol.GetSubstructMatches(ether_pattern)) if product_mol else 0
            )

            if product_ethers > reactant_ethers:
                ether_formations += 1
                print(f"Ether formation detected at depth {depth}")

        # Traverse children
        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    # Start traversal
    dfs_traverse(route)

    # Check if the strategy is present
    strategy_present = alkyl_bromide_linker_used and ether_formations >= 2

    print(f"Alkyl linker fragment coupling strategy detection:")
    print(f"  1,3-dibromopropane linker used: {alkyl_bromide_linker_used}")
    print(f"  Ether formations: {ether_formations}")
    print(f"  Strategy present: {strategy_present}")

    return strategy_present
