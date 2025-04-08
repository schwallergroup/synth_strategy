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
    Detects if the synthesis uses a methyl ester as a protecting group for a carboxylic acid,
    which is later hydrolyzed.
    """
    # Track if we found relevant transformations
    found_ester_hydrolysis = False

    def is_methyl_ester_hydrolysis(reaction_smiles):
        """Check if a reaction is a methyl ester hydrolysis"""
        # Split reaction into reactants and products
        parts = reaction_smiles.split(">")
        if len(parts) != 3:
            return False

        reactants = parts[0].split(".")
        product = parts[2]

        # Check for methyl ester in reactants
        has_methyl_ester = False
        for reactant in reactants:
            try:
                mol = Chem.MolFromSmiles(reactant)
                if mol and mol.HasSubstructMatch(Chem.MolFromSmarts("[C]-[O]-[C](=[O])-[c]")):
                    has_methyl_ester = True
                    break
            except:
                continue

        # Check for carboxylic acid in product
        try:
            prod_mol = Chem.MolFromSmiles(product)
            has_carboxylic_acid = prod_mol and prod_mol.HasSubstructMatch(
                Chem.MolFromSmarts("[O]-[C](=[O])-[c]")
            )

            # If reactant has methyl ester and product has carboxylic acid, likely an ester hydrolysis
            if has_methyl_ester and has_carboxylic_acid:
                return True
        except:
            pass

        return False

    def dfs_traverse(node, depth=0):
        nonlocal found_ester_hydrolysis

        if node["type"] == "reaction":
            # Get reaction SMILES
            if "rsmi" in node.get("metadata", {}):
                reaction_smiles = node["metadata"]["rsmi"]

                # Check if this is a methyl ester hydrolysis
                if is_methyl_ester_hydrolysis(reaction_smiles):
                    found_ester_hydrolysis = True
                    print(f"Found methyl ester hydrolysis at depth {depth}")

        # Traverse children
        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    # Start traversal
    dfs_traverse(route)

    # Check if we found methyl ester hydrolysis
    if found_ester_hydrolysis:
        print("Detected methyl ester protecting group strategy")
        return True

    return False
