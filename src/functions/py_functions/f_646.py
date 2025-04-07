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
    This function detects a Michael addition (conjugate addition) of an amine
    to an α,β-unsaturated ester as a key C-N bond forming step.
    """
    michael_addition_detected = False

    def dfs_traverse(node):
        nonlocal michael_addition_detected

        if node["type"] == "reaction":
            if "rsmi" in node.get("metadata", {}):
                rsmi = node["metadata"]["rsmi"]
                reactants = rsmi.split(">")[0].split(".")
                product = rsmi.split(">")[-1]

                # Check for α,β-unsaturated ester in reactants
                unsaturated_ester_pattern = Chem.MolFromSmarts("C=CC(=O)O[#6]")

                # Check for amine in reactants
                amine_pattern = Chem.MolFromSmarts("[NH2]")

                # Check for saturated ester with adjacent nitrogen in product
                product_pattern = Chem.MolFromSmarts("CCC(=O)O[#6]")

                # Convert SMILES to RDKit molecules
                reactant_mols = [Chem.MolFromSmiles(r) for r in reactants if r]
                product_mol = Chem.MolFromSmiles(product) if product else None

                if (
                    product_mol
                    and any(
                        mol and mol.HasSubstructMatch(unsaturated_ester_pattern)
                        for mol in reactant_mols
                    )
                    and any(
                        mol and mol.HasSubstructMatch(amine_pattern)
                        for mol in reactant_mols
                    )
                    and product_mol.HasSubstructMatch(product_pattern)
                ):
                    print("Detected Michael addition of amine to α,β-unsaturated ester")
                    michael_addition_detected = True

        # Traverse children
        for child in node.get("children", []):
            dfs_traverse(child)

    # Start traversal
    dfs_traverse(route)

    return michael_addition_detected
