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
    This function detects if the synthesis uses N-alkylation to introduce
    a functionalized side chain to a heterocyclic system.
    """
    n_alkylation_found = False

    def dfs_traverse(node, depth=0):
        nonlocal n_alkylation_found

        if node["type"] == "reaction":
            # Extract reactants and product
            rsmi = node["metadata"]["rsmi"]
            reactants = rsmi.split(">")[0].split(".")
            product = rsmi.split(">")[-1]

            # Check for N-alkylation pattern
            # Look for heterocyclic NH and alkyl halide in reactants
            nh_heterocycle_pattern = Chem.MolFromSmarts("[nH]")
            alkyl_halide_pattern = Chem.MolFromSmarts("[C][Cl,Br,I]")

            nh_present = False
            alkyl_halide_present = False

            for reactant in reactants:
                mol = Chem.MolFromSmiles(reactant)
                if mol:
                    if mol.HasSubstructMatch(nh_heterocycle_pattern):
                        nh_present = True
                    if mol.HasSubstructMatch(alkyl_halide_pattern):
                        alkyl_halide_present = True

            # Check if product has N-alkyl bond
            n_alkyl_pattern = Chem.MolFromSmarts("[n][C]")
            product_mol = Chem.MolFromSmiles(product)
            n_alkyl_present = product_mol and product_mol.HasSubstructMatch(n_alkyl_pattern)

            if (nh_present or alkyl_halide_present) and n_alkyl_present:
                print(f"N-alkylation detected at depth {depth}")
                n_alkylation_found = True

        # Traverse children
        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    # Start traversal
    dfs_traverse(route)

    return n_alkylation_found
