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
    Detects if the synthetic route employs N-alkylation of a heterocycle
    using an alkyl halide.
    """
    found_n_alkylation = False

    def dfs_traverse(node, depth=0):
        nonlocal found_n_alkylation

        if node["type"] == "reaction":
            if "metadata" in node and "rsmi" in node["metadata"]:
                rsmi = node["metadata"]["rsmi"]
                reactants = rsmi.split(">")[0].split(".")
                product = rsmi.split(">")[-1]

                # Check for N-alkylation pattern
                nh_heterocycle_pattern = Chem.MolFromSmarts("[nH]")
                alkyl_halide_pattern = Chem.MolFromSmarts("[C]-[#9,#17,#35,#53]")
                n_alkylated_pattern = Chem.MolFromSmarts("[n]-[C]")

                has_nh_heterocycle = False
                has_alkyl_halide = False

                for reactant in reactants:
                    try:
                        mol = Chem.MolFromSmiles(reactant)
                        if mol:
                            if mol.HasSubstructMatch(nh_heterocycle_pattern):
                                has_nh_heterocycle = True
                            if mol.HasSubstructMatch(alkyl_halide_pattern):
                                has_alkyl_halide = True
                    except:
                        continue

                # Check if product has N-alkylated heterocycle
                try:
                    prod_mol = Chem.MolFromSmiles(product)
                    if prod_mol and prod_mol.HasSubstructMatch(n_alkylated_pattern):
                        if has_nh_heterocycle and has_alkyl_halide:
                            found_n_alkylation = True
                            print(f"Found N-alkylation at depth {depth}")
                except:
                    pass

        # Continue traversing the tree
        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    # Start traversal from the root
    dfs_traverse(route)

    return found_n_alkylation
