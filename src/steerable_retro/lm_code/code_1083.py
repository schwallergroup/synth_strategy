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
    This function detects if the synthetic route uses N-arylation to form
    C-N bonds between aryl halides and nitrogen heterocycles.
    """
    has_n_arylation = False

    def dfs_traverse(node):
        nonlocal has_n_arylation

        if node["type"] == "reaction":
            if "metadata" in node and "rsmi" in node["metadata"]:
                rsmi = node["metadata"]["rsmi"]
                reactants = rsmi.split(">")[0].split(".")
                product = rsmi.split(">")[-1]

                # Check for aryl halide and nitrogen heterocycle in reactants
                aryl_halide_pattern = Chem.MolFromSmarts("[c]-[#9,#17,#35,#53]")
                n_heterocycle_pattern = Chem.MolFromSmarts("[#7;R]")

                has_aryl_halide = False
                has_n_heterocycle = False

                for reactant in reactants:
                    try:
                        mol = Chem.MolFromSmiles(reactant)
                        if mol and mol.HasSubstructMatch(aryl_halide_pattern):
                            has_aryl_halide = True
                        if mol and mol.HasSubstructMatch(n_heterocycle_pattern):
                            has_n_heterocycle = True
                    except:
                        continue

                # Check for aryl-N bond in product
                try:
                    product_mol = Chem.MolFromSmiles(product)
                    if product_mol:
                        aryl_n_bond_pattern = Chem.MolFromSmarts("[c]-[#7;R]")
                        if product_mol.HasSubstructMatch(aryl_n_bond_pattern):
                            if has_aryl_halide and has_n_heterocycle:
                                has_n_arylation = True
                                print(f"Detected N-arylation: {rsmi}")
                except:
                    pass

        for child in node.get("children", []):
            dfs_traverse(child)

    dfs_traverse(route)
    return has_n_arylation
