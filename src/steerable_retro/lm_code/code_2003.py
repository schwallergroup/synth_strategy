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
    This function detects if a Heck coupling (aryl halide + alkene) is present in the synthesis.
    """
    heck_coupling_found = False

    def dfs_traverse(node):
        nonlocal heck_coupling_found

        if node["type"] == "reaction":
            if "metadata" in node and "rsmi" in node["metadata"]:
                rsmi = node["metadata"]["rsmi"]
                reactants = rsmi.split(">")[0].split(".")
                product = rsmi.split(">")[-1]

                # Check for aryl halide in reactants
                aryl_halide_pattern = Chem.MolFromSmarts("[c]-[Br,Cl,I,F]")
                # Check for alkene in reactants
                alkene_pattern = Chem.MolFromSmarts("[C]=[C]")
                # Check for extended conjugation in product
                extended_conjugation = Chem.MolFromSmarts("[c]-[C]=[C]-[C](=O)")

                aryl_halide_found = False
                alkene_found = False

                for reactant in reactants:
                    reactant_mol = Chem.MolFromSmiles(reactant)
                    if not reactant_mol:
                        continue

                    if reactant_mol.HasSubstructMatch(aryl_halide_pattern):
                        aryl_halide_found = True

                    if reactant_mol.HasSubstructMatch(alkene_pattern):
                        alkene_found = True

                product_mol = Chem.MolFromSmiles(product)
                if (
                    aryl_halide_found
                    and alkene_found
                    and product_mol
                    and product_mol.HasSubstructMatch(extended_conjugation)
                ):
                    heck_coupling_found = True
                    print("Found Heck coupling (aryl halide + alkene)")

        # Traverse children
        for child in node.get("children", []):
            dfs_traverse(child)

    # Start traversal
    dfs_traverse(route)

    return heck_coupling_found
