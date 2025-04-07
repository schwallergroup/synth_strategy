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
    Detects halogen exchange or manipulation strategy in the synthesis.
    """
    found_halogen_exchange = False

    def dfs_traverse(node):
        nonlocal found_halogen_exchange

        if node["type"] == "reaction" and node.get("metadata", {}).get("rsmi"):
            rsmi = node["metadata"]["rsmi"]
            reactants_smiles = rsmi.split(">")[0].split(".")
            product_smiles = rsmi.split(">")[-1]

            try:
                reactants = [Chem.MolFromSmiles(r) for r in reactants_smiles]
                product = Chem.MolFromSmiles(product_smiles)

                if product:
                    # Check for halogen patterns
                    halogen_patterns = {
                        "Br": Chem.MolFromSmarts("c[Br]"),
                        "I": Chem.MolFromSmarts("c[I]"),
                        "Cl": Chem.MolFromSmarts("c[Cl]"),
                    }

                    # Check if product has a different halogen than reactants
                    reactant_halogens = set()
                    for mol in reactants:
                        if mol:
                            for halogen, pattern in halogen_patterns.items():
                                if mol.HasSubstructMatch(pattern):
                                    reactant_halogens.add(halogen)

                    product_halogens = set()
                    for halogen, pattern in halogen_patterns.items():
                        if product.HasSubstructMatch(pattern):
                            product_halogens.add(halogen)

                    if product_halogens and product_halogens != reactant_halogens:
                        print("Found halogen exchange or manipulation")
                        found_halogen_exchange = True
            except:
                print("Error processing reaction SMILES for halogen exchange detection")

        # Process children
        for child in node.get("children", []):
            dfs_traverse(child)

    # Start traversal
    dfs_traverse(route)
    return found_halogen_exchange
