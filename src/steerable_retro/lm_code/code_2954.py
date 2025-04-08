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
    This function detects if the synthetic route contains an aryl ether formation via SNAr with a fluoroarene.
    """
    snar_reactions = []

    def dfs_traverse(node):
        nonlocal snar_reactions

        if node["type"] == "reaction":
            if "rsmi" in node.get("metadata", {}):
                rsmi = node["metadata"]["rsmi"]
                reactants = rsmi.split(">")[0]
                product = rsmi.split(">")[-1]

                reactant_mols = [Chem.MolFromSmiles(r) for r in reactants.split(".")]
                product_mol = Chem.MolFromSmiles(product)

                # Patterns for SNAr reaction
                fluoroarene_pattern = Chem.MolFromSmarts("[F][c]")
                alcohol_pattern = Chem.MolFromSmarts("[O;H1]")
                aryl_ether_pattern = Chem.MolFromSmarts("[c][O][C]")

                has_fluoroarene = any(
                    mol and mol.HasSubstructMatch(fluoroarene_pattern) for mol in reactant_mols
                )
                has_alcohol = any(
                    mol and mol.HasSubstructMatch(alcohol_pattern) for mol in reactant_mols
                )
                has_aryl_ether = product_mol and product_mol.HasSubstructMatch(aryl_ether_pattern)

                if has_fluoroarene and has_alcohol and has_aryl_ether:
                    snar_reactions.append(node)
                    print(f"Found SNAr aryl ether formation: {rsmi}")

        # Traverse children
        for child in node.get("children", []):
            dfs_traverse(child)

    # Start traversal
    dfs_traverse(route)

    result = len(snar_reactions) > 0
    print(f"Aryl ether formation via SNAr detected: {result}")
    return result
