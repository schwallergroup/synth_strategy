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
    This function detects if the synthetic route involves N-alkylation
    with a bromoacetate fragment.
    """
    n_alkylation_found = False

    def dfs_traverse(node):
        nonlocal n_alkylation_found

        if node["type"] == "reaction":
            if "metadata" in node and "rsmi" in node["metadata"]:
                rsmi = node["metadata"]["rsmi"]
                reactants = rsmi.split(">")[0].split(".")
                product = rsmi.split(">")[-1]

                # Check for bromoacetate pattern in reactants
                bromoacetate_pattern = Chem.MolFromSmarts("Br[CH2][C](=[O])[O]")
                # Check for N-alkylated product pattern
                n_alkylated_pattern = Chem.MolFromSmarts("[#7]-[CH2][C](=[O])[O]")

                try:
                    # Check if any reactant contains bromoacetate
                    bromoacetate_present = False
                    for reactant in reactants:
                        reactant_mol = Chem.MolFromSmiles(reactant)
                        if reactant_mol and reactant_mol.HasSubstructMatch(bromoacetate_pattern):
                            bromoacetate_present = True
                            break

                    product_mol = Chem.MolFromSmiles(product)

                    if (
                        bromoacetate_present
                        and product_mol
                        and product_mol.HasSubstructMatch(n_alkylated_pattern)
                    ):
                        print(f"Found N-alkylation with bromoacetate in reaction: {rsmi}")
                        n_alkylation_found = True
                except:
                    print(f"Error processing SMILES in reaction: {rsmi}")

        for child in node.get("children", []):
            dfs_traverse(child)

    dfs_traverse(route)
    return n_alkylation_found
