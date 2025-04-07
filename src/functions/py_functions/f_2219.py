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
    This function detects a strategy involving sulfonamide formation.
    """
    has_sulfonamide_formation = False

    def dfs_traverse(node, depth=0):
        nonlocal has_sulfonamide_formation

        if node["type"] == "reaction":
            if "rsmi" in node.get("metadata", {}):
                rsmi = node["metadata"]["rsmi"]
                reactants = rsmi.split(">")[0].split(".")
                product = rsmi.split(">")[-1]

                # Check for sulfonamide formation
                has_sulfonyl_chloride = False
                has_amine = False

                for reactant in reactants:
                    react_mol = Chem.MolFromSmiles(reactant)
                    if react_mol:
                        if react_mol.HasSubstructMatch(
                            Chem.MolFromSmarts("[S](=O)(=O)[Cl]")
                        ):
                            has_sulfonyl_chloride = True
                        if react_mol.HasSubstructMatch(
                            Chem.MolFromSmarts("[N;!$(N=*)]")
                        ):
                            has_amine = True

                prod_mol = Chem.MolFromSmiles(product)
                if (
                    has_sulfonyl_chloride
                    and has_amine
                    and prod_mol
                    and prod_mol.HasSubstructMatch(Chem.MolFromSmarts("[N][S](=O)(=O)"))
                ):
                    has_sulfonamide_formation = True
                    print(f"Found sulfonamide formation at depth {depth}")

        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    dfs_traverse(route)

    if has_sulfonamide_formation:
        print("Detected sulfonamide formation strategy")

    return has_sulfonamide_formation
