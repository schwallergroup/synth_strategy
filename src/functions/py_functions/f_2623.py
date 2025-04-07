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
    This function detects a linear synthesis strategy involving addition of a PEG-like linker
    to a heterocycle.
    """
    peg_linker_addition = False

    def dfs_traverse(node):
        nonlocal peg_linker_addition

        if node["type"] == "reaction":
            # Extract reactants and product
            rsmi = node["metadata"]["rsmi"]
            reactants_smiles = rsmi.split(">")[0].split(".")
            product_smiles = rsmi.split(">")[-1]

            product_mol = Chem.MolFromSmiles(product_smiles)

            if product_mol:
                # Check for PEG-like linker pattern
                peg_pattern = Chem.MolFromSmarts("[#7]-[#6]-[#8]-[#6]-[#6]-[#8]-[#6]")
                if product_mol.HasSubstructMatch(peg_pattern):
                    # Check if product has a heterocycle
                    has_heterocycle = False
                    for ring in Chem.GetSSSR(product_mol):
                        ring_atoms = list(ring)
                        if any(
                            product_mol.GetAtomWithIdx(i).GetSymbol() != "C"
                            for i in ring_atoms
                        ):
                            has_heterocycle = True
                            break

                    if has_heterocycle:
                        peg_linker_addition = True
                        print(
                            f"Detected PEG-like linker addition to heterocycle: {rsmi}"
                        )

        # Traverse children
        for child in node.get("children", []):
            dfs_traverse(child)

    # Start traversal
    dfs_traverse(route)

    return peg_linker_addition
