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
    This function detects if the synthetic route involves a mid-stage amide bond formation.
    """
    amide_formation_found = False
    amide_formation_depth = -1
    max_depth = -1

    def dfs_traverse(node, depth=0):
        nonlocal amide_formation_found, amide_formation_depth, max_depth

        max_depth = max(max_depth, depth)

        if node["type"] == "reaction":
            # Extract reactants and product
            rsmi = node["metadata"]["rsmi"]
            reactants_smiles = rsmi.split(">")[0].split(".")
            product_smiles = rsmi.split(">")[-1]

            product_mol = Chem.MolFromSmiles(product_smiles)

            # Check for amide formation
            if product_mol:
                amide_pattern = Chem.MolFromSmarts("[#6](=[O])[#7]")
                acid_pattern = Chem.MolFromSmarts("[#6](=[O])[O]")
                amine_pattern = Chem.MolFromSmarts("[#7;!$(N=*);!$(N#*)]")

                if product_mol.HasSubstructMatch(amide_pattern):
                    # Check if reactants contain acid and amine
                    has_acid = False
                    has_amine = False

                    for r_smiles in reactants_smiles:
                        r_mol = Chem.MolFromSmiles(r_smiles)
                        if r_mol:
                            if r_mol.HasSubstructMatch(acid_pattern):
                                has_acid = True
                            if r_mol.HasSubstructMatch(amine_pattern):
                                has_amine = True

                    if has_acid and has_amine:
                        amide_formation_found = True
                        amide_formation_depth = depth
                        print(f"Amide formation found at depth {depth}")

        # Traverse children
        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    # Start traversal
    dfs_traverse(route)

    # Consider it mid-stage if it's in the middle third of the synthesis
    return amide_formation_found and max_depth / 3 < amide_formation_depth < 2 * max_depth / 3
