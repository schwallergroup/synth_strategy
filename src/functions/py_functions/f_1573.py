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
    Detects if the synthesis route involves Suzuki coupling (aryl-aryl C-C bond formation).
    """
    has_suzuki = False

    def dfs_traverse(node, depth=0):
        nonlocal has_suzuki

        if node["type"] == "reaction":
            # Get reactants and product
            rsmi = node["metadata"]["rsmi"]
            reactants_smiles = rsmi.split(">")[0].split(".")
            product_smiles = rsmi.split(">")[-1]

            # Look for boronic acid/ester pattern in reactants
            boronic_pattern = Chem.MolFromSmarts("[#6]-[#5](-[#8])-[#8]")

            # Check if any reactant has a boronic acid/ester group
            reactant_mols = [Chem.MolFromSmiles(smi) for smi in reactants_smiles]
            has_boronic = any(
                [
                    mol and mol.HasSubstructMatch(boronic_pattern)
                    for mol in reactant_mols
                ]
            )

            # Check if any reactant has a halogen (Cl, Br, I) on an aromatic ring
            halogen_pattern = Chem.MolFromSmarts("c-[#17,#35,#53]")
            has_aryl_halide = any(
                [
                    mol and mol.HasSubstructMatch(halogen_pattern)
                    for mol in reactant_mols
                ]
            )

            if has_boronic and has_aryl_halide:
                has_suzuki = True
                print(f"Detected Suzuki coupling at depth {depth}")

        # Traverse children
        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    dfs_traverse(route)
    return has_suzuki
