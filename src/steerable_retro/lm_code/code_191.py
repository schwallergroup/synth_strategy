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
    This function detects a strategy involving aromatic C-N bond formation (SNAr).
    """
    has_aromatic_c_n_formation = False

    def dfs_traverse(node):
        nonlocal has_aromatic_c_n_formation

        if node["type"] == "reaction":
            if "rsmi" in node.get("metadata", {}):
                rsmi = node["metadata"]["rsmi"]
                reactants = rsmi.split(">")[0].split(".")
                product = rsmi.split(">")[-1]

                # Check for aromatic C-N bond formation
                product_mol = Chem.MolFromSmiles(product)

                if product_mol:
                    # Look for aromatic C-N bond in product
                    aromatic_c_n_pattern = Chem.MolFromSmarts("[c][N]")

                    if product_mol.HasSubstructMatch(aromatic_c_n_pattern):
                        # Check if reactants include aryl halide and amine/ammonia
                        has_aryl_halide = False
                        has_amine = False

                        for reactant in reactants:
                            reactant_mol = Chem.MolFromSmiles(reactant)
                            if reactant_mol:
                                # Check for aryl halide
                                if reactant_mol.HasSubstructMatch(
                                    Chem.MolFromSmarts("[c][Cl,Br,I,F]")
                                ):
                                    has_aryl_halide = True
                                # Check for amine/ammonia
                                if reactant_mol.HasSubstructMatch(Chem.MolFromSmarts("[N]")):
                                    has_amine = True

                        if has_aryl_halide and has_amine:
                            has_aromatic_c_n_formation = True
                            print(f"Detected aromatic C-N bond formation: {rsmi}")

        # Traverse children
        for child in node.get("children", []):
            dfs_traverse(child)

    # Start traversal
    dfs_traverse(route)

    return has_aromatic_c_n_formation
