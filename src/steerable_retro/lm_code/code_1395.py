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
    Detects if the synthesis route involves formation of a diaryl ether (Ar-O-Ar) linkage.
    """
    has_diaryl_ether_formation = False

    def dfs_traverse(node):
        nonlocal has_diaryl_ether_formation

        if node["type"] == "reaction":
            # Extract reactants and product
            rsmi = node["metadata"]["rsmi"]
            reactants_smiles = rsmi.split(">")[0].split(".")
            product_smiles = rsmi.split(">")[-1]

            # Check for diaryl ether formation
            reactant_mols = [Chem.MolFromSmiles(smi) for smi in reactants_smiles]
            product_mol = Chem.MolFromSmiles(product_smiles)

            if product_mol and all(m for m in reactant_mols):
                # Look for diaryl ether pattern in product
                diaryl_ether_pattern = Chem.MolFromSmarts("[c][O][c]")
                if product_mol.HasSubstructMatch(diaryl_ether_pattern):
                    # Check if this pattern exists in any single reactant
                    exists_in_single_reactant = False
                    for r_mol in reactant_mols:
                        if r_mol.HasSubstructMatch(diaryl_ether_pattern):
                            exists_in_single_reactant = True
                            break

                    # If not in any single reactant, it was formed in this reaction
                    if not exists_in_single_reactant:
                        has_diaryl_ether_formation = True
                        print(f"Diaryl ether formation detected in reaction: {rsmi}")

        # Traverse children
        for child in node.get("children", []):
            dfs_traverse(child)

    # Start traversal
    dfs_traverse(route)

    print(f"Diaryl ether formation strategy: {has_diaryl_ether_formation}")
    return has_diaryl_ether_formation
