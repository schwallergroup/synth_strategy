#!/bin/python

"""LM-defined function for strategy description."""

from rdkit.Chem import AllChem, rdFMCS
import copy
from collections import deque
import rdkit.Chem as Chem
from rdkit.Chem import rdMolDescriptors
from rdkit.Chem import rdChemReactions
from rdkit.Chem import AllChem
from rdkit.Chem import rdFMCS
import rdkit.Chem.rdFMCS
from rdkit.Chem import AllChem, Descriptors, rdMolDescriptors
from rdkit import Chem
from rdkit.Chem import Descriptors
from rdkit.Chem import AllChem, rdMolDescriptors
from rdkit.Chem import AllChem, Descriptors, Lipinski
from rdkit.Chem import rdmolops
import re
from rdkit.Chem.Scaffolds import MurckoScaffold
from rdkit.Chem import AllChem, Descriptors
import traceback
import rdkit
from collections import Counter


def main(route):
    """
    Detects if the synthesis involves an N-arylation disconnection between
    an aryl group and a nitrogen heterocycle (like imidazole).
    """
    n_arylation_detected = False

    def dfs_traverse(node, depth=0):
        nonlocal n_arylation_detected

        if node["type"] == "reaction":
            rsmi = node["metadata"]["rsmi"]
            reactants = rsmi.split(">")[0].split(".")
            product = rsmi.split(">")[-1]

            product_mol = Chem.MolFromSmiles(product)

            if product_mol:
                # Check for aryl-N-heterocycle pattern in product
                aryl_n_heterocycle_pattern = Chem.MolFromSmarts("[c][N][c]1[n][c][c][c][n]1")
                if product_mol.HasSubstructMatch(aryl_n_heterocycle_pattern):
                    # Check if reactants contain aryl amine and heterocycle
                    aryl_amine_pattern = Chem.MolFromSmarts("[c][NH2]")
                    heterocycle_pattern = Chem.MolFromSmarts("[c]1[n][c][c][c][n]1")

                    has_aryl_amine = any(
                        Chem.MolFromSmiles(r)
                        and Chem.MolFromSmiles(r).HasSubstructMatch(aryl_amine_pattern)
                        for r in reactants
                        if r
                    )
                    has_heterocycle = any(
                        Chem.MolFromSmiles(r)
                        and Chem.MolFromSmiles(r).HasSubstructMatch(heterocycle_pattern)
                        for r in reactants
                        if r
                    )

                    if has_aryl_amine and has_heterocycle:
                        n_arylation_detected = True
                        print(f"N-arylation disconnection detected at depth {depth}")

        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    dfs_traverse(route)
    return n_arylation_detected
