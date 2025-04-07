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
    This function detects a synthetic strategy involving nitro group reduction to amines.
    It looks for at least one reaction where a nitro group is reduced to an amine.
    """
    nitro_to_amine_count = 0

    def dfs_traverse(node):
        nonlocal nitro_to_amine_count

        if node["type"] == "reaction" and "metadata" in node and "rsmi" in node["metadata"]:
            rsmi = node["metadata"]["rsmi"]
            reactants_smiles = rsmi.split(">")[0].split(".")
            product_smiles = rsmi.split(">")[-1]

            # Check for nitro group in reactants
            reactants_with_nitro = 0
            for reactant in reactants_smiles:
                reactant_mol = Chem.MolFromSmiles(reactant)
                if reactant_mol and reactant_mol.HasSubstructMatch(
                    Chem.MolFromSmarts("[#6]-[N+](=[O])-[O-]")
                ):
                    reactants_with_nitro += 1

            # Check for amine group in product
            product_mol = Chem.MolFromSmiles(product_smiles)
            if product_mol and product_mol.HasSubstructMatch(Chem.MolFromSmarts("[#6]-[NH2]")):
                # If reactants had nitro and product has amine, likely a reduction
                if reactants_with_nitro > 0:
                    nitro_to_amine_count += 1
                    print(f"Found nitro reduction to amine in reaction: {rsmi}")

        # Traverse children
        for child in node.get("children", []):
            dfs_traverse(child)

    # Start traversal
    dfs_traverse(route)

    # Return True if at least one nitro reduction was found
    result = nitro_to_amine_count > 0
    print(
        f"Nitro to amine interconversion strategy detected: {result} (count: {nitro_to_amine_count})"
    )
    return result
