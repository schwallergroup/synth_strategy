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
    Detects a strategy where the core scaffold (e.g., pyridine-piperazine connection)
    is assembled early in the synthesis, followed by elaboration.
    """
    pyridine_pattern = Chem.MolFromSmarts("[#6]1:[#6]:[#6]:[#7]:[#6]:[#6]:1")
    piperazine_pattern = Chem.MolFromSmarts("[#7]1[#6][#6][#7][#6][#6]1")
    early_scaffold_assembly = False

    def dfs_traverse(node, depth=0):
        nonlocal early_scaffold_assembly

        if node["type"] == "reaction" and depth >= 4:  # Early in synthesis (high depth)
            if "metadata" in node and "rsmi" in node["metadata"]:
                rsmi = node["metadata"]["rsmi"]
                reactants = rsmi.split(">")[0].split(".")
                product = rsmi.split(">")[-1]

                try:
                    prod_mol = Chem.MolFromSmiles(product)

                    # Check if product contains both pyridine and piperazine
                    if (
                        prod_mol
                        and prod_mol.HasSubstructMatch(pyridine_pattern)
                        and prod_mol.HasSubstructMatch(piperazine_pattern)
                    ):

                        # Check if reactants were separate fragments
                        pyridine_in_reactants = False
                        piperazine_in_reactants = False

                        for reactant in reactants:
                            react_mol = Chem.MolFromSmiles(reactant)
                            if react_mol:
                                if react_mol.HasSubstructMatch(pyridine_pattern):
                                    pyridine_in_reactants = True
                                if react_mol.HasSubstructMatch(piperazine_pattern):
                                    piperazine_in_reactants = True

                        # If both fragments were in reactants and now connected in product
                        if pyridine_in_reactants and piperazine_in_reactants:
                            print(f"Found early scaffold assembly at depth {depth}")
                            early_scaffold_assembly = True
                except:
                    print("Error in processing molecules for scaffold assembly")

        # Continue traversal
        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    # Start traversal from root
    dfs_traverse(route)
    return early_scaffold_assembly
