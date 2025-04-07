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
    Detects if the synthesis route uses an SNAr coupling strategy to form C-N bonds
    between an amine and a chloro-aromatic system.
    """
    snar_coupling_found = False

    def dfs_traverse(node):
        nonlocal snar_coupling_found

        if node["type"] == "reaction":
            if "rsmi" in node.get("metadata", {}):
                rsmi = node["metadata"]["rsmi"]
                reactants = rsmi.split(">")[0].split(".")
                product = rsmi.split(">")[-1]

                # Check for SNAr pattern: chloro-aromatic + amine → aryl-amine
                product_mol = Chem.MolFromSmiles(product)

                # Look for patterns in reactants
                chloro_aromatic_found = False
                amine_found = False

                for reactant in reactants:
                    reactant_mol = Chem.MolFromSmiles(reactant)
                    if not reactant_mol:
                        continue

                    # Check for chloro-aromatic
                    chloro_aromatic_pattern = Chem.MolFromSmarts("[Cl][c]")
                    if reactant_mol.HasSubstructMatch(chloro_aromatic_pattern):
                        chloro_aromatic_found = True

                    # Check for amine
                    amine_pattern = Chem.MolFromSmarts("[N;!$(N=*);!$(NC=O)]")
                    if reactant_mol.HasSubstructMatch(amine_pattern):
                        amine_found = True

                # Check if product has aryl-amine bond
                if product_mol:
                    aryl_amine_pattern = Chem.MolFromSmarts("[c][N;!$(N=*);!$(NC=O)]")
                    if (
                        chloro_aromatic_found
                        and amine_found
                        and product_mol.HasSubstructMatch(aryl_amine_pattern)
                    ):
                        snar_coupling_found = True
                        print("SNAr coupling detected")

        # Traverse children
        for child in node.get("children", []):
            dfs_traverse(child)

    dfs_traverse(route)
    return snar_coupling_found
