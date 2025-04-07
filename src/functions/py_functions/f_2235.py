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
    Detects if the synthesis uses a reductive amination strategy
    (aldehyde/ketone + amine → amine)
    """
    has_reductive_amination = False

    def dfs_traverse(node, depth=0):
        nonlocal has_reductive_amination

        if node["type"] == "reaction":
            if "rsmi" in node.get("metadata", {}):
                rsmi = node["metadata"]["rsmi"]
                reactants = rsmi.split(">")[0].split(".")
                product = rsmi.split(">")[-1]

                # Check for reductive amination pattern
                has_aldehyde_or_ketone = False
                has_amine = False

                for reactant in reactants:
                    r_mol = Chem.MolFromSmiles(reactant)
                    if r_mol:
                        # Check for aldehyde or ketone
                        if r_mol.HasSubstructMatch(
                            Chem.MolFromSmarts("[#6](=[#8])[#6,#1]")
                        ):
                            has_aldehyde_or_ketone = True
                        # Check for amine
                        if r_mol.HasSubstructMatch(
                            Chem.MolFromSmarts("[NH,NH2;!$(NC=O)]")
                        ):
                            has_amine = True

                # Check if product has a new C-N bond
                if has_aldehyde_or_ketone and has_amine:
                    prod_mol = Chem.MolFromSmiles(product)
                    if prod_mol and prod_mol.HasSubstructMatch(
                        Chem.MolFromSmarts("[#6]-[#7;!$(NC=O)]")
                    ):
                        print(f"Found reductive amination at depth {depth}")
                        has_reductive_amination = True

        # Traverse children
        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    # Start traversal from root
    dfs_traverse(route)
    return has_reductive_amination
