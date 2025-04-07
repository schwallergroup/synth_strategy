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
    This function detects if the synthesis route involves the construction of a complex
    nitrogen-heterocyclic system.
    """
    heterocyclic_construction = False

    def dfs_traverse(node):
        nonlocal heterocyclic_construction

        if (
            node["type"] == "reaction"
            and "metadata" in node
            and "rsmi" in node["metadata"]
        ):
            rsmi = node["metadata"]["rsmi"]
            reactants = rsmi.split(">")[0].split(".")
            product = rsmi.split(">")[-1]

            # Check if product contains multiple nitrogen-containing rings
            product_mol = Chem.MolFromSmiles(product)

            if product_mol:
                # Count nitrogen atoms in rings
                n_in_rings = 0
                for atom in product_mol.GetAtoms():
                    if atom.GetAtomicNum() == 7 and atom.IsInRing():
                        n_in_rings += 1

                # If product has multiple nitrogens in rings, check if they were formed in this step
                if n_in_rings >= 2:
                    # Check if any reactant has fewer nitrogen-containing rings
                    for reactant in reactants:
                        reactant_mol = Chem.MolFromSmiles(reactant)
                        if reactant_mol:
                            reactant_n_in_rings = 0
                            for atom in reactant_mol.GetAtoms():
                                if atom.GetAtomicNum() == 7 and atom.IsInRing():
                                    reactant_n_in_rings += 1

                            if reactant_n_in_rings < n_in_rings:
                                print(f"Heterocyclic construction detected: {rsmi}")
                                heterocyclic_construction = True
                                break

        # Continue traversing
        for child in node.get("children", []):
            dfs_traverse(child)

    # Start traversal from the root
    dfs_traverse(route)
    return heterocyclic_construction
