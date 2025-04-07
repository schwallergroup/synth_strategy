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
    This function detects late-stage alcohol mesylation (conversion of alcohol to mesylate
    in the final steps of the synthesis)
    """
    # Track if mesylation was found and at what depth
    mesylation_found = False
    mesylation_depth = float("inf")

    def dfs_traverse(node, depth=0):
        nonlocal mesylation_found, mesylation_depth

        if node["type"] == "reaction":
            # Extract reactants and product
            rsmi = node["metadata"]["rsmi"]
            reactants_smiles = rsmi.split(">")[0].split(".")
            product_smiles = rsmi.split(">")[-1]

            reactant_mols = [Chem.MolFromSmiles(r) for r in reactants_smiles]
            product_mol = Chem.MolFromSmiles(product_smiles)

            # Check for alcohol in reactants
            alcohol_pattern = Chem.MolFromSmarts("[#6][#8;H1]")

            # Check for mesylate in product
            mesylate_pattern = Chem.MolFromSmarts("[#6][#8][#16](=[#8])(=[#8])[#6]")

            # Check for mesylation: alcohol to mesylate
            if any(
                mol is not None and mol.HasSubstructMatch(alcohol_pattern)
                for mol in reactant_mols
            ) and product_mol.HasSubstructMatch(mesylate_pattern):
                mesylation_found = True
                mesylation_depth = min(mesylation_depth, depth)

        # Traverse children
        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    # Start traversal
    dfs_traverse(route)

    # Check if mesylation was found in the late stage (depth <= 1)
    late_stage = mesylation_found and mesylation_depth <= 1

    print(f"Late-stage alcohol mesylation detected: {late_stage}")
    print(f"Mesylation depth: {mesylation_depth if mesylation_found else 'Not found'}")

    return late_stage
