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
    This function detects if the synthetic route contains an aryl amine formation reaction.
    """
    aryl_amine_reactions = []

    def dfs_traverse(node):
        nonlocal aryl_amine_reactions

        if node["type"] == "reaction":
            if "rsmi" in node.get("metadata", {}):
                rsmi = node["metadata"]["rsmi"]
                reactants = rsmi.split(">")[0]
                product = rsmi.split(">")[-1]

                reactant_mols = [Chem.MolFromSmiles(r) for r in reactants.split(".")]
                product_mol = Chem.MolFromSmiles(product)

                # Patterns for aryl amine formation
                aryl_halide_pattern = Chem.MolFromSmarts("[c][Br,I,Cl]")
                amine_pattern = Chem.MolFromSmarts("[N;H2]")
                aryl_amine_pattern = Chem.MolFromSmarts("[c][N;!$(N=*);!$(N#*)]")

                has_aryl_halide = any(
                    mol and mol.HasSubstructMatch(aryl_halide_pattern)
                    for mol in reactant_mols
                )
                has_amine = any(
                    mol and mol.HasSubstructMatch(amine_pattern)
                    for mol in reactant_mols
                )
                has_aryl_amine = product_mol and product_mol.HasSubstructMatch(
                    aryl_amine_pattern
                )

                if has_aryl_halide and has_amine and has_aryl_amine:
                    aryl_amine_reactions.append(node)
                    print(f"Found aryl amine formation reaction: {rsmi}")

        # Traverse children
        for child in node.get("children", []):
            dfs_traverse(child)

    # Start traversal
    dfs_traverse(route)

    result = len(aryl_amine_reactions) > 0
    print(f"Aryl amine formation detected: {result}")
    return result
