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
    Detects if the synthesis uses N-alkylation with benzyl halides.
    """
    has_benzyl_alkylation = False

    def dfs_traverse(node, depth=0):
        nonlocal has_benzyl_alkylation

        if node["type"] == "reaction":
            if "rsmi" in node.get("metadata", {}):
                rsmi = node["metadata"]["rsmi"]
                reactants = rsmi.split(">")[0].split(".")
                product = rsmi.split(">")[-1]

                # Check for benzyl halide pattern in reactants
                benzyl_halide_pattern = Chem.MolFromSmarts("c1ccccc1C[Br,Cl,I,F]")
                amine_pattern = Chem.MolFromSmarts("[NH]")

                has_benzyl_halide = False
                has_amine = False

                for reactant in reactants:
                    try:
                        mol = Chem.MolFromSmiles(reactant)
                        if mol:
                            if mol.HasSubstructMatch(benzyl_halide_pattern):
                                has_benzyl_halide = True
                            if mol.HasSubstructMatch(amine_pattern):
                                has_amine = True
                    except:
                        continue

                # Check if product has new C-N bond with benzyl group
                if has_benzyl_halide and has_amine:
                    try:
                        product_mol = Chem.MolFromSmiles(product)
                        if product_mol:
                            benzyl_amine_pattern = Chem.MolFromSmarts("c1ccccc1C[N]")
                            if product_mol.HasSubstructMatch(benzyl_amine_pattern):
                                has_benzyl_alkylation = True
                                print(f"Found benzyl alkylation at depth {depth}")
                    except:
                        pass

        # Traverse children
        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    dfs_traverse(route)
    return has_benzyl_alkylation
