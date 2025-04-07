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
    Detects if the synthetic route employs a Boc protection strategy
    for an amine group.
    """
    # Track if we found Boc protection
    found_boc_protection = False

    def is_boc_protection(reaction_smiles):
        """Check if a reaction is a Boc protection of an amine"""
        # Split into reactants and product
        parts = reaction_smiles.split(">")
        if len(parts) < 3:
            return False

        reactants = parts[0].split(".")
        product = parts[2]

        # Check for amine in reactants
        amine_pattern = Chem.MolFromSmarts("[#7;H2]")

        # Check for Boc group in product
        boc_pattern = Chem.MolFromSmarts("CC(C)(C)OC(=O)[#7]")

        has_amine = False
        for reactant in reactants:
            try:
                mol = Chem.MolFromSmiles(reactant)
                if mol and mol.HasSubstructMatch(amine_pattern):
                    has_amine = True
                    break
            except:
                continue

        has_boc = False
        try:
            prod_mol = Chem.MolFromSmiles(product)
            if prod_mol and prod_mol.HasSubstructMatch(boc_pattern):
                has_boc = True
        except:
            pass

        return has_amine and has_boc

    def dfs_traverse(node):
        nonlocal found_boc_protection

        # Check if this is a reaction node
        if node.get("type") == "reaction":
            # Get reaction SMILES from metadata
            if "rsmi" in node.get("metadata", {}):
                rsmi = node["metadata"]["rsmi"]

                # Check if this is a Boc protection
                if is_boc_protection(rsmi):
                    found_boc_protection = True
                    print("Found Boc protection reaction")

        # Traverse children
        for child in node.get("children", []):
            dfs_traverse(child)

    # Start traversal
    dfs_traverse(route)

    return found_boc_protection
