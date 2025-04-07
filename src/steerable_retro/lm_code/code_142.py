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
    Detects if the synthesis route includes a Williamson ether synthesis
    (formation of ether from phenol and alkyl halide)
    """
    found_williamson = False

    def dfs_traverse(node):
        nonlocal found_williamson

        if node["type"] == "reaction":
            if "rsmi" in node.get("metadata", {}):
                rsmi = node["metadata"]["rsmi"]
                reactants = rsmi.split(">")[0].split(".")
                product = rsmi.split(">")[-1]

                # Create RDKit mol objects
                product_mol = Chem.MolFromSmiles(product)
                reactant_mols = [Chem.MolFromSmiles(r) for r in reactants if r]

                if product_mol and all(reactant_mols):
                    # SMARTS for phenol
                    phenol_pattern = Chem.MolFromSmarts("[c][OH]")

                    # SMARTS for alkyl halide
                    alkyl_halide_pattern = Chem.MolFromSmarts("[#6][CH2][Br,Cl,I]")

                    # SMARTS for aryl-alkyl ether
                    ether_pattern = Chem.MolFromSmarts("[c][O][CH2][#6]")

                    # Check if reactants contain phenol and alkyl halide
                    has_phenol = any(r.HasSubstructMatch(phenol_pattern) for r in reactant_mols)
                    has_alkyl_halide = any(
                        r.HasSubstructMatch(alkyl_halide_pattern) for r in reactant_mols
                    )

                    # Check if product has ether
                    has_ether = product_mol.HasSubstructMatch(ether_pattern)

                    if has_phenol and has_alkyl_halide and has_ether:
                        print("Found Williamson ether synthesis")
                        found_williamson = True

        # Continue traversal
        for child in node.get("children", []):
            dfs_traverse(child)

    # Start traversal from the root
    dfs_traverse(route)
    return found_williamson
