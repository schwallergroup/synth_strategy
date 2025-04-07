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
    This function detects if the synthetic route involves reduction of an alkene (C=C to C-C).
    """
    alkene_pattern = Chem.MolFromSmarts("C=C")
    has_alkene_reduction = False

    def dfs_traverse(node):
        nonlocal has_alkene_reduction

        if node["type"] == "reaction":
            if "metadata" in node and "rsmi" in node["metadata"]:
                rsmi = node["metadata"]["rsmi"]
                reactants_smiles = rsmi.split(">")[0]
                product_smiles = rsmi.split(">")[-1]

                try:
                    reactants_mol = Chem.MolFromSmiles(reactants_smiles)
                    product_mol = Chem.MolFromSmiles(product_smiles)

                    # Check if reactants have alkene but product doesn't
                    if (
                        reactants_mol
                        and reactants_mol.HasSubstructMatch(alkene_pattern)
                        and (not product_mol or not product_mol.HasSubstructMatch(alkene_pattern))
                    ):
                        print(f"Detected alkene reduction in reaction: {rsmi}")
                        has_alkene_reduction = True
                except:
                    print(f"Error processing reaction SMILES: {rsmi}")

        for child in node.get("children", []):
            dfs_traverse(child)

    dfs_traverse(route)
    return has_alkene_reduction
