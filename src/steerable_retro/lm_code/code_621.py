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
    This function detects if the synthetic route involves multiple aromatic nucleophilic
    substitution reactions where a halogen is replaced by O or N.
    """
    substitution_count = 0

    def dfs_traverse(node):
        nonlocal substitution_count

        if node["type"] == "reaction":
            if "rsmi" in node.get("metadata", {}):
                rsmi = node["metadata"]["rsmi"]
                reactants_smiles = rsmi.split(">")[0]
                product_smiles = rsmi.split(">")[-1]

                try:
                    reactants = [Chem.MolFromSmiles(r) for r in reactants_smiles.split(".")]
                    product = Chem.MolFromSmiles(product_smiles)

                    if all(r is not None for r in reactants) and product is not None:
                        # Check for aryl halide in reactants
                        aryl_halide_pattern = Chem.MolFromSmarts("[c]-[F,Cl,Br,I]")
                        # Check for aryl ether or aryl amine in product
                        aryl_ether_pattern = Chem.MolFromSmarts("[c]-[O]-[#6]")
                        aryl_amine_pattern = Chem.MolFromSmarts("[c]-[N]")

                        reactants_have_aryl_halide = any(
                            r.HasSubstructMatch(aryl_halide_pattern) for r in reactants
                        )
                        product_has_aryl_ether = product.HasSubstructMatch(aryl_ether_pattern)
                        product_has_aryl_amine = product.HasSubstructMatch(aryl_amine_pattern)

                        if reactants_have_aryl_halide and (
                            product_has_aryl_ether or product_has_aryl_amine
                        ):
                            print(f"Aromatic substitution detected in reaction: {rsmi}")
                            substitution_count += 1
                except:
                    print(f"Error processing reaction SMILES: {rsmi}")

        # Process children
        for child in node.get("children", []):
            dfs_traverse(child)

    # Start traversal from the root
    dfs_traverse(route)

    # Return True if multiple aromatic substitutions are detected
    return substitution_count >= 2
