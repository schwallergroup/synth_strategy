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
    This function detects a synthetic strategy involving the construction of a bicyclic heterocycle
    system from a monocyclic precursor.
    """
    bicyclic_construction_detected = False

    def dfs_traverse(node):
        nonlocal bicyclic_construction_detected

        if node["type"] == "reaction":
            if "metadata" in node and "rsmi" in node["metadata"]:
                rsmi = node["metadata"]["rsmi"]
                reactants_smiles = rsmi.split(">")[0].split(".")
                product_smiles = rsmi.split(">")[-1]

                # Convert to RDKit molecules
                reactants_mols = [Chem.MolFromSmiles(r) for r in reactants_smiles if r]
                product_mol = Chem.MolFromSmiles(product_smiles)

                if all(reactants_mols) and product_mol:
                    # Check for monocyclic heterocycle in reactants
                    monocyclic_pattern = Chem.MolFromSmarts(
                        "[#6]1[#6,#7,#8,#16][#6,#7,#8,#16][#6,#7,#8,#16][#6,#7,#8,#16][#6,#7,#8,#16]1"
                    )

                    # Check for bicyclic heterocycle in product
                    # This is a simplified pattern - in practice you'd need more specific patterns
                    bicyclic_pattern = Chem.MolFromSmarts(
                        "[#6]1[#6,#7,#8,#16][#6,#7,#8,#16][#6]2[#6,#7,#8,#16][#6,#7,#8,#16][#6,#7,#8,#16][#6,#7,#8,#16]2[#6,#7,#8,#16]1"
                    )

                    reactants_monocyclic = any(
                        [
                            mol.HasSubstructMatch(monocyclic_pattern)
                            for mol in reactants_mols
                        ]
                    )
                    product_bicyclic = product_mol.HasSubstructMatch(bicyclic_pattern)

                    if reactants_monocyclic and product_bicyclic:
                        print("Bicyclic heterocycle construction detected")
                        bicyclic_construction_detected = True

        # Traverse children
        for child in node.get("children", []):
            dfs_traverse(child)

    # Start traversal from the root
    dfs_traverse(route)
    return bicyclic_construction_detected
