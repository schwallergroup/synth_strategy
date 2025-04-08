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
    Detects chalcone cleavage as a key disconnection strategy.
    Looks for breaking of α,β-unsaturated ketone into aldehyde and ketone.
    """
    chalcone_cleavage_found = False

    def dfs_traverse(node, depth=0):
        nonlocal chalcone_cleavage_found

        if node["type"] == "reaction":
            if "rsmi" in node.get("metadata", {}):
                rsmi = node["metadata"]["rsmi"]
                reactants = rsmi.split(">")[0].split(".")
                product = rsmi.split(">")[-1]

                # Check if product has chalcone structure (C=C-C=O connected to aromatics)
                product_mol = Chem.MolFromSmiles(product)
                if product_mol:
                    chalcone_pattern = Chem.MolFromSmarts("c-/C=C/C(=O)c")
                    if product_mol.HasSubstructMatch(chalcone_pattern):
                        # Check if reactants contain aldehyde and ketone
                        has_aldehyde = False
                        has_ketone = False

                        for reactant in reactants:
                            r_mol = Chem.MolFromSmiles(reactant)
                            if r_mol:
                                aldehyde_pattern = Chem.MolFromSmarts("c-C=O")
                                if r_mol.HasSubstructMatch(aldehyde_pattern):
                                    has_aldehyde = True

                                ketone_pattern = Chem.MolFromSmarts("C-C(=O)-c")
                                if r_mol.HasSubstructMatch(ketone_pattern):
                                    has_ketone = True

                        if has_aldehyde and has_ketone:
                            print("Chalcone cleavage detected at depth", depth)
                            chalcone_cleavage_found = True

        # Traverse children
        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    # Start traversal
    dfs_traverse(route)

    return chalcone_cleavage_found
