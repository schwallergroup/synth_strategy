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
    Detects if the synthesis involves isoxazole formation via cycloaddition.
    """
    result = False

    def dfs_traverse(node):
        nonlocal result

        if (
            node["type"] == "reaction"
            and "metadata" in node
            and "rsmi" in node["metadata"]
        ):
            rsmi = node["metadata"]["rsmi"]
            reactants = rsmi.split(">")[0].split(".")
            product = rsmi.split(">")[-1]

            try:
                product_mol = Chem.MolFromSmiles(product)
                reactant_mols = [Chem.MolFromSmiles(r) for r in reactants]

                # Check for isoxazole pattern in product
                isoxazole_pattern = Chem.MolFromSmarts("[c]1[n][o][c][c]1")
                if product_mol and product_mol.HasSubstructMatch(isoxazole_pattern):
                    # Check if isoxazole is newly formed
                    isoxazole_exists_in_reactants = False
                    for r_mol in reactant_mols:
                        if r_mol and r_mol.HasSubstructMatch(isoxazole_pattern):
                            isoxazole_exists_in_reactants = True
                            break

                    if not isoxazole_exists_in_reactants:
                        # Check if one reactant has an alkyne and another has an oxime
                        has_alkyne = False
                        has_oxime = False

                        alkyne_pattern = Chem.MolFromSmarts("[C]#[C]")
                        oxime_pattern = Chem.MolFromSmarts("[C]=[N][OH]")

                        for r_mol in reactant_mols:
                            if r_mol:
                                if r_mol.HasSubstructMatch(alkyne_pattern):
                                    has_alkyne = True
                                if r_mol.HasSubstructMatch(oxime_pattern):
                                    has_oxime = True

                        if has_alkyne and has_oxime:
                            print("Detected isoxazole formation via cycloaddition")
                            result = True
            except:
                pass

        # Continue traversing
        for child in node.get("children", []):
            dfs_traverse(child)

    # Start traversal
    dfs_traverse(route)
    return result
