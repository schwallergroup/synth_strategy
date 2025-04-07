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
    Detects if the synthesis route involves formation of heterocyclic rings,
    particularly pyrimidine formation.
    """
    heterocycle_formed = False

    def dfs_traverse(node):
        nonlocal heterocycle_formed

        if node["type"] == "reaction":
            if "rsmi" in node["metadata"]:
                rsmi = node["metadata"]["rsmi"]
                reactants = rsmi.split(">")[0].split(".")
                product = rsmi.split(">")[-1]

                try:
                    # Count rings in reactants and product
                    reactant_mols = [
                        Chem.MolFromSmiles(r) for r in reactants if Chem.MolFromSmiles(r)
                    ]
                    product_mol = Chem.MolFromSmiles(product)

                    if product_mol and reactant_mols:
                        # Count heterocyclic rings in reactants
                        reactant_rings = 0
                        for mol in reactant_mols:
                            if mol:
                                # Count rings containing N, O, or S
                                hetero_pattern = Chem.MolFromSmarts(
                                    "[#7,#8,#16]~1~[#6]~[#6]~[#6]~[#6,#7]~1"
                                )
                                reactant_rings += len(mol.GetSubstructMatches(hetero_pattern))

                        # Count heterocyclic rings in product
                        product_rings = 0
                        if product_mol:
                            hetero_pattern = Chem.MolFromSmarts(
                                "[#7,#8,#16]~1~[#6]~[#6]~[#6]~[#6,#7]~1"
                            )
                            product_rings = len(product_mol.GetSubstructMatches(hetero_pattern))

                            # Specifically check for pyrimidine formation
                            pyrimidine_pattern = Chem.MolFromSmarts("[#7]1[#6][#7][#6][#6][#6]1")
                            if product_mol.HasSubstructMatch(pyrimidine_pattern):
                                print(f"Pyrimidine ring detected in product: {product}")
                                heterocycle_formed = True

                        # Check if new heterocyclic rings were formed
                        if product_rings > reactant_rings:
                            heterocycle_formed = True
                            print(
                                f"Heterocycle formation detected: {reactant_rings} -> {product_rings} rings"
                            )
                except:
                    pass

        # Continue traversing
        for child in node.get("children", []):
            dfs_traverse(child)

    # Start traversal
    dfs_traverse(route)

    print(f"Heterocycle formation strategy detected: {heterocycle_formed}")
    return heterocycle_formed
