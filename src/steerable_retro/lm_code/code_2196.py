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
    Detects a sequence of sulfur oxidation state changes (thiol→thioether→sulfone→thioether)
    in the synthetic route.
    """
    # Track sulfur oxidation states through the route
    thiol_found = False
    thioether_found = False
    sulfone_found = False
    sulfone_to_thioether = False

    def dfs_traverse(node):
        nonlocal thiol_found, thioether_found, sulfone_found, sulfone_to_thioether

        if node["type"] == "reaction" and "metadata" in node and "rsmi" in node["metadata"]:
            rsmi = node["metadata"]["rsmi"]
            reactants = rsmi.split(">")[0].split(".")
            product = rsmi.split(">")[-1]

            # Check for thiol in reactants
            thiol_pattern = Chem.MolFromSmarts("[SH]")
            for reactant in reactants:
                try:
                    mol = Chem.MolFromSmiles(reactant)
                    if mol and mol.HasSubstructMatch(thiol_pattern):
                        thiol_found = True
                        print(f"Found thiol in reactant: {reactant}")
                except:
                    continue

            # Check for thioether in products
            thioether_pattern = Chem.MolFromSmarts("[#16]-[#6]")
            try:
                mol = Chem.MolFromSmiles(product)
                if mol and mol.HasSubstructMatch(thioether_pattern):
                    thioether_found = True
                    print(f"Found thioether in product: {product}")
            except:
                pass

            # Check for sulfone in products
            sulfone_pattern = Chem.MolFromSmarts("[#16](=[#8])(=[#8])-[#6]")
            try:
                mol = Chem.MolFromSmiles(product)
                if mol and mol.HasSubstructMatch(sulfone_pattern):
                    sulfone_found = True
                    print(f"Found sulfone in product: {product}")
            except:
                pass

            # Check for sulfone to thioether transformation
            if any(
                Chem.MolFromSmiles(r) and Chem.MolFromSmiles(r).HasSubstructMatch(sulfone_pattern)
                for r in reactants
                if Chem.MolFromSmiles(r)
            ):
                try:
                    mol = Chem.MolFromSmiles(product)
                    if (
                        mol
                        and mol.HasSubstructMatch(thioether_pattern)
                        and not mol.HasSubstructMatch(sulfone_pattern)
                    ):
                        sulfone_to_thioether = True
                        print(f"Found sulfone to thioether transformation")
                except:
                    pass

        # Process children
        for child in node.get("children", []):
            dfs_traverse(child)

    # Start traversal
    dfs_traverse(route)

    # Return True if we found the sequence thiol→thioether→sulfone→thioether
    result = thiol_found and thioether_found and sulfone_found and sulfone_to_thioether
    print(f"Sulfur oxidation state changes detected: {result}")
    return result
