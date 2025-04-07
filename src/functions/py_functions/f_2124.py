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
    This function detects a synthetic strategy involving formation of aromatic ethers.
    It looks for reactions that form C-O bonds between aromatic rings and other groups.
    """
    ether_formations = 0

    def dfs_traverse(node):
        nonlocal ether_formations

        if (
            node["type"] == "reaction"
            and "metadata" in node
            and "rsmi" in node["metadata"]
        ):
            rsmi = node["metadata"]["rsmi"]
            reactants_smiles = rsmi.split(">")[0].split(".")
            product_smiles = rsmi.split(">")[-1]

            # Check for aromatic ether in product
            product_mol = Chem.MolFromSmiles(product_smiles)
            if product_mol and product_mol.HasSubstructMatch(
                Chem.MolFromSmarts("c-[O]-c")
            ):
                # Check if this bond was formed in this reaction (not present in reactants)
                ether_in_reactants = False
                for reactant in reactants_smiles:
                    reactant_mol = Chem.MolFromSmiles(reactant)
                    if reactant_mol and reactant_mol.HasSubstructMatch(
                        Chem.MolFromSmarts("c-[O]-c")
                    ):
                        ether_in_reactants = True
                        break

                if not ether_in_reactants:
                    ether_formations += 1
                    print(f"Found aromatic ether formation: {rsmi}")

        # Traverse children
        for child in node.get("children", []):
            dfs_traverse(child)

    # Start traversal
    dfs_traverse(route)

    # Return True if at least one aromatic ether formation was found
    result = ether_formations > 0
    print(
        f"Aromatic ether formation strategy detected: {result} (count: {ether_formations})"
    )
    return result
