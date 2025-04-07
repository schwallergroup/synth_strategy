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
    This function detects if the synthetic route involves multiple amide bond formations.
    """
    amide_formation_count = 0

    def dfs_traverse(node):
        nonlocal amide_formation_count

        if (
            node["type"] == "reaction"
            and "metadata" in node
            and "rsmi" in node["metadata"]
        ):
            rsmi = node["metadata"]["rsmi"]
            reactants_smiles = rsmi.split(">")[0]
            product_smiles = rsmi.split(">")[-1]

            try:
                reactants_mol = Chem.MolFromSmiles(reactants_smiles)
                product_mol = Chem.MolFromSmiles(product_smiles)

                if reactants_mol and product_mol:
                    # Check for amide formation
                    amide_pattern = Chem.MolFromSmarts("[NX3][CX3](=[OX1])")

                    reactants_amide_count = len(
                        reactants_mol.GetSubstructMatches(amide_pattern)
                    )
                    product_amide_count = len(
                        product_mol.GetSubstructMatches(amide_pattern)
                    )

                    if product_amide_count > reactants_amide_count:
                        print(f"Amide formation detected: {rsmi}")
                        amide_formation_count += 1
            except:
                print(f"Error processing reaction SMILES: {rsmi}")

        for child in node.get("children", []):
            dfs_traverse(child)

    dfs_traverse(route)
    return (
        amide_formation_count >= 2
    )  # Return True if at least 2 amide formations detected
