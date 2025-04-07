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
    This function detects if the synthetic route involves interconversion
    between aldehydes and vinyl ethers.
    """
    # Track if we found the pattern
    found_interconversion = False

    def dfs_traverse(node):
        nonlocal found_interconversion

        if (
            node["type"] == "reaction"
            and "metadata" in node
            and "rsmi" in node["metadata"]
        ):
            rsmi = node["metadata"]["rsmi"]
            reactants = rsmi.split(">")[0].split(".")
            product = rsmi.split(">")[-1]

            # Check for aldehyde pattern
            aldehyde_pattern = "[#6]-[#6H1]=[O]"
            # Check for vinyl ether pattern
            vinyl_ether_pattern = "[#6]-[#8]-[#6]=[#6]"

            # Check for aldehyde to vinyl ether conversion
            reactant_has_aldehyde = False
            for reactant in reactants:
                mol = Chem.MolFromSmiles(reactant)
                if mol and mol.HasSubstructMatch(Chem.MolFromSmarts(aldehyde_pattern)):
                    reactant_has_aldehyde = True
                    break

            product_mol = Chem.MolFromSmiles(product)
            product_has_vinyl_ether = product_mol and product_mol.HasSubstructMatch(
                Chem.MolFromSmarts(vinyl_ether_pattern)
            )

            if reactant_has_aldehyde and product_has_vinyl_ether:
                print("Found aldehyde to vinyl ether conversion")
                found_interconversion = True

            # Check for vinyl ether to aldehyde conversion
            reactant_has_vinyl_ether = False
            for reactant in reactants:
                mol = Chem.MolFromSmiles(reactant)
                if mol and mol.HasSubstructMatch(
                    Chem.MolFromSmarts(vinyl_ether_pattern)
                ):
                    reactant_has_vinyl_ether = True
                    break

            product_has_aldehyde = product_mol and product_mol.HasSubstructMatch(
                Chem.MolFromSmarts(aldehyde_pattern)
            )

            if reactant_has_vinyl_ether and product_has_aldehyde:
                print("Found vinyl ether to aldehyde conversion")
                found_interconversion = True

        # Continue traversing
        for child in node.get("children", []):
            dfs_traverse(child)

    # Start traversal
    dfs_traverse(route)
    return found_interconversion
