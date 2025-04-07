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
    This function detects if the synthetic route contains multiple C-N bond formations.
    """
    cn_bond_count = 0

    def dfs_traverse(node):
        nonlocal cn_bond_count

        if node["type"] == "reaction":
            if "rsmi" in node.get("metadata", {}):
                rsmi = node["metadata"]["rsmi"]
                reactants_smiles = rsmi.split(">")[0]
                product_smiles = rsmi.split(">")[-1]

                # Check for amide formation
                reactants = Chem.MolFromSmiles(reactants_smiles)
                product = Chem.MolFromSmiles(product_smiles)

                if reactants and product:
                    # Look for carboxylic acid or derivative in reactants
                    acid_pattern = Chem.MolFromSmarts("[C](=[O])[O,N]")
                    amine_pattern = Chem.MolFromSmarts("[N;!$(N=*);!$(NC=O)]")
                    amide_pattern = Chem.MolFromSmarts("[N;!$(N=*)][C](=[O])")

                    has_acid = reactants.HasSubstructMatch(acid_pattern)
                    has_amine = reactants.HasSubstructMatch(amine_pattern)
                    has_amide_product = product.HasSubstructMatch(amide_pattern)

                    # Check for C-N bond formation (including amide)
                    if has_acid and has_amine and has_amide_product:
                        cn_bond_count += 1
                        print(f"C-N bond formation detected in reaction: {rsmi}")

                    # Check for other C-N bond formations (e.g., SNAr)
                    aryl_halide = Chem.MolFromSmarts("[c][Br,Cl,I,F]")
                    if reactants.HasSubstructMatch(aryl_halide) and reactants.HasSubstructMatch(
                        amine_pattern
                    ):
                        # Check if product has new C-N bond where halide was
                        if product.HasSubstructMatch(Chem.MolFromSmarts("[c][N]")):
                            cn_bond_count += 1
                            print(f"C-N bond formation (SNAr) detected in reaction: {rsmi}")

        for child in node.get("children", []):
            dfs_traverse(child)

    dfs_traverse(route)
    return cn_bond_count >= 2
