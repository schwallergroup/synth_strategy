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
    Detects if the synthesis route involves reduction of a nitro group to an amine.
    """
    nitro_reduction_detected = False

    def dfs_traverse(node, depth=0):
        nonlocal nitro_reduction_detected

        if node["type"] == "reaction":
            rsmi = node["metadata"]["rsmi"]
            reactants_smiles = rsmi.split(">")[0].split(".")
            product_smiles = rsmi.split(">")[-1]

            # Create RDKit molecules
            try:
                reactant_mols = [Chem.MolFromSmiles(r) for r in reactants_smiles]
                product_mol = Chem.MolFromSmiles(product_smiles)

                if product_mol and all(r for r in reactant_mols):
                    # Check for nitro group in reactants
                    nitro_pattern = Chem.MolFromSmarts("[N+](=O)[O-]")
                    amine_pattern = Chem.MolFromSmarts("[NH2]")

                    nitro_in_reactants = False
                    for r_mol in reactant_mols:
                        if r_mol.HasSubstructMatch(nitro_pattern):
                            nitro_in_reactants = True
                            break

                    # Check for amine in product
                    amine_in_product = product_mol.HasSubstructMatch(amine_pattern)

                    if nitro_in_reactants and amine_in_product:
                        nitro_reduction_detected = True
                        print(f"Nitro to amine reduction detected at depth {depth}")
                        return
            except:
                print(f"Error processing reaction at depth {depth}")

        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    dfs_traverse(route)
    return nitro_reduction_detected
