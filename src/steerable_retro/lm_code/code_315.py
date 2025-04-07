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
    This function detects a synthetic strategy involving N-arylation
    (formation of C-N bond between aryl halide and amine).
    """
    n_arylation_found = False

    def dfs_traverse(node, depth=0):
        nonlocal n_arylation_found

        if node["type"] == "reaction":
            if "rsmi" in node.get("metadata", {}):
                rsmi = node["metadata"]["rsmi"]
                reactants = rsmi.split(">")[0].split(".")
                product = rsmi.split(">")[-1]

                # Check for secondary amine in reactants
                amine_pattern = Chem.MolFromSmarts("[#7;H1]")

                # Check for aryl halide in reactants
                aryl_halide_pattern = Chem.MolFromSmarts("[c]-[Br,I,Cl]")

                # Check for aryl-N bond in product that wasn't in reactants
                aryl_n_pattern = Chem.MolFromSmarts("[c]-[#7]")

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

                # Check if product has a new C-N bond
                if has_amine and has_aryl_halide:
                    try:
                        product_mol = Chem.MolFromSmiles(product)
                        if product_mol and product_mol.HasSubstructMatch(aryl_n_pattern):
                            print(f"N-arylation detected at depth {depth}")
                            n_arylation_found = True
                    except:
                        pass

        # Traverse children
        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    # Start traversal from the root
    dfs_traverse(route)

    return n_arylation_found
