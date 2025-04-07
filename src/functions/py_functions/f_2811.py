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
    Detects if the synthetic route involves sequential nitrogen functionalization steps
    (multiple reactions forming C-N bonds)
    """
    n_functionalization_depths = []

    def dfs_traverse(node, depth=0):
        if (
            node["type"] == "reaction"
            and "metadata" in node
            and "rsmi" in node["metadata"]
        ):
            rsmi = node["metadata"]["rsmi"]
            reactants_part = rsmi.split(">")[0]
            product_part = rsmi.split(">")[-1]

            # Check for nitrogen-containing functional groups
            amine_pattern = Chem.MolFromSmarts("[N;H]")
            amide_pattern = Chem.MolFromSmarts("C(=O)[N]")

            # Check if product has a new C-N bond not present in reactants
            try:
                product_mol = Chem.MolFromSmiles(product_part)

                # Look for amide formation or amine alkylation
                if product_mol and (
                    product_mol.HasSubstructMatch(amide_pattern)
                    or product_mol.HasSubstructMatch(
                        Chem.MolFromSmarts("[#6]-[N](-[#6])-[#6]")
                    )
                ):

                    # Check if this is a new bond formation
                    reactants = reactants_part.split(".")
                    has_same_bond_in_reactants = False

                    for reactant in reactants:
                        try:
                            r_mol = Chem.MolFromSmiles(reactant)
                            if r_mol and r_mol.HasSubstructMatch(amide_pattern):
                                has_same_bond_in_reactants = True
                        except:
                            continue

                    if not has_same_bond_in_reactants:
                        print(f"Found nitrogen functionalization at depth {depth}")
                        n_functionalization_depths.append(depth)
            except:
                pass

        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    dfs_traverse(route)

    # Check if we have at least 2 nitrogen functionalization steps
    if len(n_functionalization_depths) >= 2:
        # Check if they are sequential (no more than 2 steps apart)
        n_functionalization_depths.sort()
        for i in range(1, len(n_functionalization_depths)):
            if n_functionalization_depths[i] - n_functionalization_depths[i - 1] <= 2:
                print("Found sequential nitrogen functionalization steps")
                return True

    return False
