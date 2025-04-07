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
    This function detects a convergent peptide synthesis strategy where multiple fragments
    are combined through amide bond formations.
    """
    convergent_steps = 0
    amide_formations = 0

    def dfs_traverse(node, depth=0):
        nonlocal convergent_steps, amide_formations

        if node["type"] == "reaction":
            if "rsmi" in node.get("metadata", {}):
                rsmi = node["metadata"]["rsmi"]
                reactants = rsmi.split(">")[0].split(".")
                product = rsmi.split(">")[-1]

                # Check if this is a convergent step (multiple reactants -> one product)
                if len(reactants) >= 2:
                    product_mol = Chem.MolFromSmiles(product)
                    reactant_mols = [Chem.MolFromSmiles(r) for r in reactants]

                    # Check if product has an amide bond
                    if product_mol and product_mol.HasSubstructMatch(
                        Chem.MolFromSmarts("[C](=[O])[NH]")
                    ):
                        # Check if this is forming a new amide bond
                        amide_count_product = len(
                            product_mol.GetSubstructMatches(
                                Chem.MolFromSmarts("[C](=[O])[NH]")
                            )
                        )
                        amide_count_reactants = sum(
                            len(
                                mol.GetSubstructMatches(
                                    Chem.MolFromSmarts("[C](=[O])[NH]")
                                )
                            )
                            for mol in reactant_mols
                            if mol
                        )

                        if amide_count_product > amide_count_reactants:
                            convergent_steps += 1
                            amide_formations += 1
                            print(f"Found convergent amide formation at depth {depth}")

        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    dfs_traverse(route)

    result = convergent_steps >= 1 and amide_formations >= 2
    print(f"Convergent peptide synthesis: {result}")
    print(f"Convergent steps: {convergent_steps}")
    print(f"Amide formations: {amide_formations}")

    return result
