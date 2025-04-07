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
    Detects if the synthesis involves an aromatic nucleophilic substitution
    to form a C-N bond between an amine and an aromatic system.
    """
    has_aromatic_amine_coupling = False

    def dfs_traverse(node, depth=0):
        nonlocal has_aromatic_amine_coupling

        if node["type"] == "reaction":
            rsmi = node["metadata"]["rsmi"]
            reactants_smiles = rsmi.split(">")[0].split(".")
            product_smiles = rsmi.split(">")[-1]

            # Look for fluoroaromatic in reactants
            has_fluoroaromatic = False
            has_amine = False

            for reactant in reactants_smiles:
                mol = Chem.MolFromSmiles(reactant)
                if mol:
                    fluoro_aromatic_pattern = Chem.MolFromSmarts("[F][c]")
                    amine_pattern = Chem.MolFromSmarts("[NX3;H2,H1;!$(NC=O)]")

                    if mol.HasSubstructMatch(fluoro_aromatic_pattern):
                        has_fluoroaromatic = True
                    if mol.HasSubstructMatch(amine_pattern):
                        has_amine = True

            # Check if product has new C-N bond to aromatic
            if has_fluoroaromatic and has_amine:
                product_mol = Chem.MolFromSmiles(product_smiles)
                if product_mol:
                    aryl_amine_pattern = Chem.MolFromSmarts("[NX3;!$(NC=O)][c]")
                    if product_mol.HasSubstructMatch(aryl_amine_pattern):
                        has_aromatic_amine_coupling = True
                        print(f"Aromatic amine coupling detected at depth {depth}")

        # Traverse children
        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    dfs_traverse(route)
    return has_aromatic_amine_coupling
