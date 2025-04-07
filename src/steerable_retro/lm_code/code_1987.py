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
    Detects if the synthesis route involves a convergent approach where complex fragments are joined.
    """
    convergent_synthesis = False

    def dfs_traverse(node, depth=0):
        nonlocal convergent_synthesis

        if node["type"] == "reaction":
            if "rsmi" in node.get("metadata", {}):
                rsmi = node["metadata"]["rsmi"]
                reactants = rsmi.split(">")[0].split(".")
                product = rsmi.split(">")[-1]

                # Count complex fragments (defined as having more than 10 atoms)
                complex_fragment_count = 0
                for reactant in reactants:
                    reactant_mol = Chem.MolFromSmiles(reactant)
                    if reactant_mol and reactant_mol.GetNumAtoms() > 10:
                        complex_fragment_count += 1

                # If we have at least 2 complex fragments being joined, it's convergent
                if complex_fragment_count >= 2:
                    product_mol = Chem.MolFromSmiles(product)
                    if product_mol:
                        # Check if product is more complex than individual reactants
                        if product_mol.GetNumAtoms() > max(
                            [
                                Chem.MolFromSmiles(r).GetNumAtoms()
                                for r in reactants
                                if Chem.MolFromSmiles(r)
                            ]
                        ):
                            print(
                                f"Convergent synthesis detected at depth {depth} with {complex_fragment_count} complex fragments"
                            )
                            convergent_synthesis = True

        # Traverse children
        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    dfs_traverse(route)
    return convergent_synthesis
