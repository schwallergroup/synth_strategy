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
    This function detects if the synthesis employs a convergent multi-component approach
    to form a heterocycle, specifically looking for reactions with 3+ reactants forming a heterocycle.
    """
    found_strategy = False

    def dfs_traverse(node, depth=0):
        nonlocal found_strategy

        if node["type"] == "reaction" and "rsmi" in node.get("metadata", {}):
            rsmi = node["metadata"]["rsmi"]
            reactants_smiles = rsmi.split(">")[0].split(".")
            product_smiles = rsmi.split(">")[-1]

            # Check if we have 3 or more reactants
            if len(reactants_smiles) >= 3:
                product_mol = Chem.MolFromSmiles(product_smiles)
                if product_mol:
                    # Check if product contains a heterocycle
                    heterocycle_pattern = Chem.MolFromSmarts(
                        "[*;!#6;!#1]~1~[*]~[*]~[*]~[*]~[*]~1"
                    )  # 6-membered heterocycle
                    pyrimidine_pattern = Chem.MolFromSmarts(
                        "[n;r6]1[c;r6][n;r6][c;r6][c;r6][c;r6]1"
                    )

                    if product_mol.HasSubstructMatch(
                        heterocycle_pattern
                    ) or product_mol.HasSubstructMatch(pyrimidine_pattern):
                        # Check if the heterocycle is not present in any of the reactants
                        reactant_has_heterocycle = False
                        for reactant in reactants_smiles:
                            reactant_mol = Chem.MolFromSmiles(reactant)
                            if reactant_mol and (
                                reactant_mol.HasSubstructMatch(heterocycle_pattern)
                                or reactant_mol.HasSubstructMatch(pyrimidine_pattern)
                            ):
                                reactant_has_heterocycle = True
                                break

                        if not reactant_has_heterocycle:
                            found_strategy = True
                            print(
                                f"Found convergent multi-component heterocycle formation at depth {depth}"
                            )

        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    dfs_traverse(route)
    return found_strategy
