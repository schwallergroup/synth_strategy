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
    Detects a strategy involving reduction of a nitro group to an amine.
    """
    found_nitro_reduction = False

    def dfs_traverse(node, depth=0):
        nonlocal found_nitro_reduction

        if node["type"] == "reaction":
            # Extract reactants and product
            rsmi = node["metadata"]["rsmi"]
            reactants_smiles = rsmi.split(">")[0].split(".")
            product_smiles = rsmi.split(">")[-1]

            # Convert to RDKit molecules
            reactant_mols = [Chem.MolFromSmiles(r) for r in reactants_smiles if r]
            product_mol = Chem.MolFromSmiles(product_smiles)

            if product_mol and all(r for r in reactant_mols):
                # Check for nitro group in reactants
                nitro_patt = Chem.MolFromSmarts("[N+](=O)[O-]")
                # Check for amine in product
                amine_patt = Chem.MolFromSmarts("[NH2]")

                reactant_has_nitro = False
                for r in reactant_mols:
                    if r.HasSubstructMatch(nitro_patt):
                        reactant_has_nitro = True
                        break

                product_has_amine = product_mol.HasSubstructMatch(amine_patt)

                # If reactant has nitro and product has amine, it's likely a nitro reduction
                if reactant_has_nitro and product_has_amine:
                    found_nitro_reduction = True
                    print(f"Found nitro reduction at depth {depth}")

        # Traverse children
        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    # Start traversal
    dfs_traverse(route)

    return found_nitro_reduction
