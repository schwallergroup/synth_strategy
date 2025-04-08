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
    This function detects if the synthetic route involves an amide formation reaction
    between a carboxylic acid and an amine.
    """
    amide_formation_detected = False

    def dfs_traverse(node):
        nonlocal amide_formation_detected

        if node["type"] == "reaction":
            if "rsmi" in node.get("metadata", {}):
                rsmi = node["metadata"]["rsmi"]
                reactants_smiles = rsmi.split(">")[0]
                product_smiles = rsmi.split(">")[-1]

                try:
                    reactants = [Chem.MolFromSmiles(r) for r in reactants_smiles.split(".")]
                    product = Chem.MolFromSmiles(product_smiles)

                    if all(r is not None for r in reactants) and product is not None:
                        # Check for carboxylic acid in reactants
                        carboxylic_acid_pattern = Chem.MolFromSmarts("[CX3](=[OX1])[OX2H1]")
                        # Check for amine in reactants
                        amine_pattern = Chem.MolFromSmarts("[NX3;H2,H1;!$(NC=O)]")
                        # Check for amide in product
                        amide_pattern = Chem.MolFromSmarts("[CX3](=[OX1])[NX3]")

                        reactants_have_acid = any(
                            r.HasSubstructMatch(carboxylic_acid_pattern) for r in reactants
                        )
                        reactants_have_amine = any(
                            r.HasSubstructMatch(amine_pattern) for r in reactants
                        )
                        product_has_amide = product.HasSubstructMatch(amide_pattern)

                        if reactants_have_acid and reactants_have_amine and product_has_amide:
                            print(f"Amide formation detected in reaction: {rsmi}")
                            amide_formation_detected = True
                except:
                    print(f"Error processing reaction SMILES: {rsmi}")

        # Process children
        for child in node.get("children", []):
            dfs_traverse(child)

    # Start traversal from the root
    dfs_traverse(route)

    return amide_formation_detected
