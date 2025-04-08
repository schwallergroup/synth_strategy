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
    This function detects a synthetic strategy where two aromatic fragments
    are connected via ether bonds using an alkyl chain linker.
    """
    # Track if we found the pattern
    found_pattern = False
    # Track the presence of key features
    has_ether_linkage = False
    has_alkyl_chain = False
    aromatic_fragments = 0

    def dfs_traverse(node):
        nonlocal found_pattern, has_ether_linkage, has_alkyl_chain, aromatic_fragments

        if node["type"] == "mol":
            # Check final product for the desired pattern
            mol = Chem.MolFromSmiles(node["smiles"])
            if mol:
                # Check for aromatic-O-alkyl-O-aromatic pattern
                aromatic_ether_pattern = Chem.MolFromSmarts(
                    "c-[#8]-[#6]-[#6]-[#6]-[#6]-[#6]-[#8]-c"
                )
                if mol.HasSubstructMatch(aromatic_ether_pattern):
                    has_ether_linkage = True
                    has_alkyl_chain = True
                    aromatic_fragments = 2
                    found_pattern = True
                    print(
                        f"Found aromatic-O-alkyl-O-aromatic pattern in molecule: {node['smiles']}"
                    )

        elif node["type"] == "reaction":
            # Extract reactants and product
            rsmi = node["metadata"].get("rsmi", "")
            if not rsmi:
                return

            reactants = rsmi.split(">")[0].split(".")
            product = rsmi.split(">")[-1]

            # Check for ether formation reactions
            if any("O" in r for r in reactants) and "O" in product:
                prod_mol = Chem.MolFromSmiles(product)
                if prod_mol and prod_mol.HasSubstructMatch(Chem.MolFromSmarts("c-[#8]-[#6]")):
                    has_ether_linkage = True
                    print(f"Found ether formation reaction: {rsmi}")

                # Check for alkyl chain
                if any(
                    Chem.MolFromSmiles(r)
                    and Chem.MolFromSmiles(r).HasSubstructMatch(
                        Chem.MolFromSmarts("[#6]-[#6]-[#6]-[#6]-[#6]")
                    )
                    for r in reactants
                    if Chem.MolFromSmiles(r)
                ):
                    has_alkyl_chain = True
                    print(f"Found alkyl chain in reaction: {rsmi}")

                # Count aromatic fragments
                for r in reactants:
                    r_mol = Chem.MolFromSmiles(r)
                    if r_mol and r_mol.HasSubstructMatch(Chem.MolFromSmarts("c1ccccc1")):
                        aromatic_fragments += 1
                        print(f"Found aromatic fragment: {r}")

        # Continue traversing
        for child in node.get("children", []):
            dfs_traverse(child)

    # Start traversal
    dfs_traverse(route)

    # Check if we have all the required elements
    if has_ether_linkage and has_alkyl_chain and aromatic_fragments >= 2:
        found_pattern = True

    return found_pattern
