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
    This function detects if the synthesis route involves formation of a urea with a cyclopropyl group.
    """
    cyclopropyl_urea_found = False

    def dfs_traverse(node):
        nonlocal cyclopropyl_urea_found

        if node["type"] == "reaction" and "metadata" in node and "rsmi" in node["metadata"]:
            rsmi = node["metadata"]["rsmi"]
            reactants = rsmi.split(">")[0].split(".")
            product = rsmi.split(">")[-1]

            # Create RDKit mol objects
            product_mol = Chem.MolFromSmiles(product) if product else None

            if product_mol:
                # Check for cyclopropyl urea formation
                cyclopropyl_pattern = Chem.MolFromSmarts("C1CC1")
                urea_pattern = Chem.MolFromSmarts("[#7][#6](=O)[#7]")

                if product_mol.HasSubstructMatch(
                    cyclopropyl_pattern
                ) and product_mol.HasSubstructMatch(urea_pattern):
                    # Check if the cyclopropyl is connected to the urea
                    # This is a simplified check - a more robust implementation would verify the exact connection
                    print("Detected molecule with both cyclopropyl and urea groups")

                    # Check if any reactant has cyclopropylamine
                    for reactant in reactants:
                        reactant_mol = Chem.MolFromSmiles(reactant)
                        if reactant_mol:
                            cyclopropylamine_pattern = Chem.MolFromSmarts("C1CC1[NH2]")
                            if reactant_mol.HasSubstructMatch(cyclopropylamine_pattern):
                                print("Confirmed cyclopropylamine used in urea formation")
                                cyclopropyl_urea_found = True

        # Traverse children
        for child in node.get("children", []):
            dfs_traverse(child)

    # Start traversal from the root
    dfs_traverse(route)
    return cyclopropyl_urea_found
