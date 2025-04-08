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
    This function specifically detects imidazole ring formation in the synthetic route.
    """
    imidazole_formation_detected = False

    def dfs_traverse(node):
        nonlocal imidazole_formation_detected

        if node["type"] == "reaction":
            rsmi = node["metadata"]["rsmi"]
            reactants_smiles = rsmi.split(">")[0].split(".")
            product_smiles = rsmi.split(">")[-1]

            # Check for imidazole pattern in product
            imidazole_pattern = Chem.MolFromSmarts("c1ncnc1")

            # Check if product has an imidazole ring
            product_mol = Chem.MolFromSmiles(product_smiles)
            if product_mol and product_mol.HasSubstructMatch(imidazole_pattern):
                # Check if reactants don't have imidazole
                imidazole_in_reactants = False
                for reactant_smiles in reactants_smiles:
                    reactant_mol = Chem.MolFromSmiles(reactant_smiles)
                    if reactant_mol and reactant_mol.HasSubstructMatch(imidazole_pattern):
                        imidazole_in_reactants = True
                        break

                if not imidazole_in_reactants:
                    imidazole_formation_detected = True
                    print(f"Detected imidazole ring formation: {rsmi}")

        for child in node.get("children", []):
            dfs_traverse(child)

    dfs_traverse(route)

    print(f"Imidazole ring formation detected: {imidazole_formation_detected}")
    return imidazole_formation_detected
