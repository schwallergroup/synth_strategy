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
    This function detects late-stage reductive amination (final or penultimate step).
    Specifically, it looks for C=O + HN → CH2-N transformation at low depth (0-1).
    """
    reductive_amination_found = False
    depth_of_reductive_amination = None

    def dfs_traverse(node, depth=0):
        nonlocal reductive_amination_found, depth_of_reductive_amination

        if node["type"] == "reaction":
            if "rsmi" in node.get("metadata", {}):
                rsmi = node["metadata"]["rsmi"]
                reactants = rsmi.split(">")[0].split(".")
                product = rsmi.split(">")[-1]

                # Check for reductive amination pattern
                # Look for aldehyde/ketone + amine in reactants and C-N bond in product
                has_carbonyl = False
                has_amine = False

                for reactant in reactants:
                    mol = Chem.MolFromSmiles(reactant)
                    if mol:
                        # Check for aldehyde or ketone
                        if mol.HasSubstructMatch(Chem.MolFromSmarts("[C;$(C=O)]")):
                            has_carbonyl = True
                        # Check for amine
                        if mol.HasSubstructMatch(
                            Chem.MolFromSmarts("[#7;!$(N=*);!$(NC=O);!$(NS=O)]")
                        ):
                            has_amine = True

                # Check if product has new C-N bond where C is not carbonyl
                if has_carbonyl and has_amine:
                    prod_mol = Chem.MolFromSmiles(product)
                    if prod_mol and prod_mol.HasSubstructMatch(
                        Chem.MolFromSmarts("[C;!$(C=O)]-[#7]")
                    ):
                        reductive_amination_found = True
                        depth_of_reductive_amination = depth
                        print(f"Reductive amination found at depth {depth}")

        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    dfs_traverse(route)

    # Check if reductive amination occurs at late stage (depth 0-1)
    if (
        reductive_amination_found
        and depth_of_reductive_amination is not None
        and depth_of_reductive_amination <= 1
    ):
        print("Late-stage reductive amination strategy detected")
        return True
    return False
