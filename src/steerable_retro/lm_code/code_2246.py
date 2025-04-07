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
    This function detects if the synthetic route involves formation of a nitrogen heterocycle.
    """
    heterocycle_formation_found = False

    def dfs_traverse(node):
        nonlocal heterocycle_formation_found

        if node["type"] == "reaction" and "metadata" in node and "rsmi" in node["metadata"]:
            rsmi = node["metadata"]["rsmi"]
            reactants = rsmi.split(">")[0].split(".")
            product = rsmi.split(">")[-1]

            # Check for nitrogen heterocycle formation
            # Patterns for common N-heterocycles
            pyrimidine_pattern = Chem.MolFromSmarts("c1ncncc1")
            pyridine_pattern = Chem.MolFromSmarts("c1ncccc1")
            imidazole_pattern = Chem.MolFromSmarts("c1ncnc1")

            product_mol = Chem.MolFromSmiles(product)

            if product_mol:
                has_heterocycle = (
                    product_mol.HasSubstructMatch(pyrimidine_pattern)
                    or product_mol.HasSubstructMatch(pyridine_pattern)
                    or product_mol.HasSubstructMatch(imidazole_pattern)
                )

                if has_heterocycle:
                    # Check if reactants don't have the same heterocycle
                    heterocycle_in_reactants = False
                    for reactant in reactants:
                        reactant_mol = Chem.MolFromSmiles(reactant)
                        if reactant_mol and (
                            reactant_mol.HasSubstructMatch(pyrimidine_pattern)
                            or reactant_mol.HasSubstructMatch(pyridine_pattern)
                            or reactant_mol.HasSubstructMatch(imidazole_pattern)
                        ):
                            heterocycle_in_reactants = True
                            break

                    if not heterocycle_in_reactants:
                        print("Nitrogen heterocycle formation detected")
                        heterocycle_formation_found = True

        # Traverse children
        for child in node.get("children", []):
            dfs_traverse(child)

    # Start traversal from root
    dfs_traverse(route)
    return heterocycle_formation_found
