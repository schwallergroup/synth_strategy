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
    This function detects a linear synthesis with sequential functionalization
    of an aromatic core structure.
    """
    is_linear = True
    aromatic_functionalization_count = 0

    def dfs_traverse(node):
        nonlocal is_linear, aromatic_functionalization_count

        if node["type"] == "reaction":
            # Extract reactants and product
            rsmi = node["metadata"]["rsmi"]
            reactants = rsmi.split(">")[0].split(".")
            product = rsmi.split(">")[-1]

            # Check if this is a linear step (1-2 reactants)
            if len(reactants) > 2:
                is_linear = False

            # Check if this involves aromatic functionalization
            try:
                product_mol = Chem.MolFromSmiles(product)
                for reactant in reactants:
                    reactant_mol = Chem.MolFromSmiles(reactant)
                    if reactant_mol and product_mol:
                        # Check if both have aromatic rings
                        if any(atom.GetIsAromatic() for atom in reactant_mol.GetAtoms()) and any(
                            atom.GetIsAromatic() for atom in product_mol.GetAtoms()
                        ):
                            # Check if product has more functional groups than reactant
                            if sum(
                                1 for atom in product_mol.GetAtoms() if atom.GetAtomicNum() > 1
                            ) > sum(
                                1 for atom in reactant_mol.GetAtoms() if atom.GetAtomicNum() > 1
                            ):
                                aromatic_functionalization_count += 1
                                break
            except:
                print("Error in processing molecules for aromatic functionalization check")

        # Process children
        for child in node.get("children", []):
            dfs_traverse(child)

    # Call dfs_traverse on the root node
    dfs_traverse(route)

    # Return True if it's a linear synthesis with at least 3 aromatic functionalization steps
    result = is_linear and aromatic_functionalization_count >= 3
    print(f"Linear synthesis with sequential aromatic functionalization: {result}")
    print(f"  - Linear synthesis: {is_linear}")
    print(f"  - Aromatic functionalization steps: {aromatic_functionalization_count}")

    return result
