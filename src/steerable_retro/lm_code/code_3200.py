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
    Detects a strategy using nucleophilic aromatic substitution for C-N bond formation.
    Looks for amine and aryl halide reactants forming a C-N bond.
    """
    snar_found = False

    def dfs_traverse(node):
        nonlocal snar_found

        if node["type"] == "reaction":
            if "rsmi" in node.get("metadata", {}):
                rsmi = node["metadata"]["rsmi"]
                reactants = rsmi.split(">")[0].split(".")
                product = rsmi.split(">")[-1]

                # Check if we have multiple reactants
                if len(reactants) >= 2:
                    # Look for amine and aryl halide patterns
                    amine_pattern = Chem.MolFromSmarts("[#7;H1,H2]")
                    aryl_halide_pattern = Chem.MolFromSmarts("[#6]:[#6]-[#35,#53,#17]")

                    has_amine = False
                    has_aryl_halide = False

                    for reactant in reactants:
                        try:
                            mol = Chem.MolFromSmiles(reactant)
                            if mol and mol.HasSubstructMatch(amine_pattern):
                                has_amine = True
                            if mol and mol.HasSubstructMatch(aryl_halide_pattern):
                                has_aryl_halide = True
                        except:
                            continue

                    # Check if both patterns are found
                    if has_amine and has_aryl_halide:
                        try:
                            prod_mol = Chem.MolFromSmiles(product)
                            # Check if product has C-N bond
                            if prod_mol and prod_mol.HasSubstructMatch(
                                Chem.MolFromSmarts("[#6]:[#6]-[#7]")
                            ):
                                print("Found potential nucleophilic aromatic substitution")
                                snar_found = True
                        except:
                            pass

        # Continue traversing
        for child in node.get("children", []):
            dfs_traverse(child)

    # Start traversal
    dfs_traverse(route)
    return snar_found
