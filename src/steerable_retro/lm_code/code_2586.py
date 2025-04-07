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
    Detects if the synthesis involves sulfonamide formation from
    a sulfonyl chloride and an amine.
    """
    sulfonamide_formation_found = False

    def dfs_traverse(node, depth=0):
        nonlocal sulfonamide_formation_found

        if node["type"] == "reaction":
            rsmi = node["metadata"]["rsmi"]
            reactants_smiles = rsmi.split(">")[0].split(".")
            product_smiles = rsmi.split(">")[-1]

            # Check for sulfonamide in product
            sulfonamide_pattern = Chem.MolFromSmarts("[#16](=[#8])(=[#8])[#7]")
            product_mol = Chem.MolFromSmiles(product_smiles)

            if product_mol and product_mol.HasSubstructMatch(sulfonamide_pattern):
                # Check for sulfonyl chloride in reactants
                sulfonyl_chloride_pattern = Chem.MolFromSmarts("[#16](=[#8])(=[#8])[Cl]")
                sulfonyl_chloride_found = False
                amine_found = False

                for reactant_smiles in reactants_smiles:
                    reactant_mol = Chem.MolFromSmiles(reactant_smiles)
                    if not reactant_mol:
                        continue

                    if reactant_mol.HasSubstructMatch(sulfonyl_chloride_pattern):
                        sulfonyl_chloride_found = True

                    # Check for amine
                    amine_pattern = Chem.MolFromSmarts("[#7;H2]")
                    if reactant_mol.HasSubstructMatch(amine_pattern):
                        amine_found = True

                if sulfonyl_chloride_found and amine_found:
                    sulfonamide_formation_found = True
                    print(f"Found sulfonamide formation at depth {depth}")

        # Traverse children
        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    # Start traversal
    dfs_traverse(route)

    if sulfonamide_formation_found:
        print("Detected sulfonamide formation strategy")
        return True
    return False
