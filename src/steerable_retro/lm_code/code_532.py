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
    Detects a synthetic strategy involving multiple nucleophilic aromatic
    substitution (SNAr) reactions.
    """
    snar_count = 0

    def dfs_traverse(node):
        nonlocal snar_count

        if node["type"] == "reaction":
            if "metadata" in node and "rsmi" in node["metadata"]:
                rsmi = node["metadata"]["rsmi"]
                reactants = rsmi.split(">")[0].split(".")
                product = rsmi.split(">")[-1]

                # Check for SNAr (nucleophilic aromatic substitution)
                aryl_halide_pattern = Chem.MolFromSmarts("[#6]:[#6]-[F,Cl,Br,I]")
                amine_pattern = Chem.MolFromSmarts("[N;H0,H1,H2]")
                aryl_amine_pattern = Chem.MolFromSmarts("[#6]:[#6]-[N;H0,H1,H2]")

                has_aryl_halide = False
                has_amine = False
                for reactant in reactants:
                    mol = Chem.MolFromSmiles(reactant)
                    if mol:
                        if mol.HasSubstructMatch(aryl_halide_pattern):
                            has_aryl_halide = True
                        if mol.HasSubstructMatch(amine_pattern):
                            has_amine = True

                product_mol = Chem.MolFromSmiles(product)
                if (
                    has_aryl_halide
                    and has_amine
                    and product_mol
                    and product_mol.HasSubstructMatch(aryl_amine_pattern)
                ):
                    snar_count += 1
                    print(f"Found SNAr reaction (count: {snar_count})")

        # Traverse children
        for child in node.get("children", []):
            dfs_traverse(child)

    # Start traversal
    dfs_traverse(route)

    # Check if multiple SNAr reactions are present
    return snar_count >= 2
