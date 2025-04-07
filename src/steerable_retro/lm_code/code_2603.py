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
    This function detects a strategy involving derivatization of a carboxylic acid group,
    particularly esterification.
    """
    has_esterification = False

    def dfs_traverse(node, depth=0):
        nonlocal has_esterification

        if node["type"] == "reaction":
            # Extract reactants and product
            rsmi = node["metadata"]["rsmi"]
            reactants_smiles = rsmi.split(">")[0].split(".")
            product_smiles = rsmi.split(">")[-1]

            # Convert to RDKit molecules
            reactants = [Chem.MolFromSmiles(r) for r in reactants_smiles if r]
            product = Chem.MolFromSmiles(product_smiles)

            if product and all(r for r in reactants):
                # Check for esterification
                acid_pattern = Chem.MolFromSmarts("[#6](=[#8])([#8;H1])")
                ester_pattern = Chem.MolFromSmarts("[#6](=[#8])([#8][#6&D1])")

                product_acid_count = len(product.GetSubstructMatches(acid_pattern))
                reactants_acid_count = sum(
                    len(r.GetSubstructMatches(acid_pattern)) for r in reactants
                )
                product_ester_count = len(product.GetSubstructMatches(ester_pattern))
                reactants_ester_count = sum(
                    len(r.GetSubstructMatches(ester_pattern)) for r in reactants
                )

                if (
                    product_ester_count > reactants_ester_count
                    and reactants_acid_count > product_acid_count
                ):
                    has_esterification = True
                    print(f"Detected esterification at depth {depth}")

                    # Check if methanol is one of the reactants
                    methanol_pattern = Chem.MolFromSmarts("[#6&D1][#8]")
                    if any(
                        r.GetSubstructMatch(methanol_pattern)
                        and Chem.MolToSmiles(r) in ["CO", "OC"]
                        for r in reactants
                    ):
                        print("Methanol used as esterification agent")

        # Traverse children
        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    # Start traversal
    dfs_traverse(route)

    print(f"Carboxylic acid derivatization strategy detected: {has_esterification}")
    return has_esterification
