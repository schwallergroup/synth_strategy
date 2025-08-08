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
    Detects N-benzylation of a nitrogen heterocycle (like tetrahydroisoquinoline)
    """
    found_n_benzylation = False

    def dfs_traverse(node):
        nonlocal found_n_benzylation

        if node["type"] == "reaction":
            rsmi = node["metadata"]["rsmi"]
            reactants = rsmi.split(">")[0].split(".")
            product = rsmi.split(">")[-1]

            # Check for N-benzylation pattern
            if len(reactants) == 2:  # Typically two reactants for N-alkylation
                product_mol = Chem.MolFromSmiles(product)

                if product_mol:
                    # SMARTS for N-benzylated heterocycle
                    n_benzyl_pattern = Chem.MolFromSmarts("[#7;R]-[#6]-[c]1[c][c][c][c][c]1")

                    # Check if the product has this pattern
                    if product_mol.HasSubstructMatch(n_benzyl_pattern):
                        # Check if one reactant has the heterocycle and one has the benzyl group
                        heterocycle_found = False
                        benzyl_found = False

                        for reactant in reactants:
                            reactant_mol = Chem.MolFromSmiles(reactant)
                            if reactant_mol:
                                # SMARTS for nitrogen heterocycle
                                heterocycle_pattern = Chem.MolFromSmarts("[#7;R]")
                                # SMARTS for benzyl halide
                                benzyl_halide_pattern = Chem.MolFromSmarts(
                                    "[#6]-[c]1[c][c][c][c][c]1"
                                )

                                if reactant_mol.HasSubstructMatch(heterocycle_pattern):
                                    heterocycle_found = True
                                if reactant_mol.HasSubstructMatch(benzyl_halide_pattern):
                                    benzyl_found = True

                        if heterocycle_found and benzyl_found:
                            found_n_benzylation = True
                            print("Found N-benzylation of nitrogen heterocycle")

        # Traverse children
        for child in node.get("children", []):
            dfs_traverse(child)

    # Start traversal
    dfs_traverse(route)

    return found_n_benzylation
