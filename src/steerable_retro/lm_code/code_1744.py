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
    This function detects if the synthetic route employs a strategy of using
    amide couplings to connect major fragments of the molecule.
    """
    # Track amide couplings
    amide_couplings = []

    def dfs_traverse(node, depth=0):
        nonlocal amide_couplings

        if node["type"] == "reaction":
            # Extract reactants and product
            rsmi = node["metadata"]["rsmi"]
            reactants_smiles = rsmi.split(">")[0].split(".")
            product_smiles = rsmi.split(">")[-1]

            # Convert to RDKit molecules
            reactants = [Chem.MolFromSmiles(smi) for smi in reactants_smiles if smi]
            product = Chem.MolFromSmiles(product_smiles) if product_smiles else None

            if product and len(reactants) >= 2:  # Need at least two reactants for fragment coupling
                # Check for carboxylic acid or activated ester in reactants
                acid_pattern = Chem.MolFromSmarts("[#6]-[C](=[O])[OH]")
                ester_pattern = Chem.MolFromSmarts("[#6]-[C](=[O])[O][#6]")
                acid_reactant = any(mol.HasSubstructMatch(acid_pattern) for mol in reactants if mol)
                ester_reactant = any(
                    mol.HasSubstructMatch(ester_pattern) for mol in reactants if mol
                )

                # Check for amine in reactants
                amine_pattern = Chem.MolFromSmarts("[#6]-[NH][#6]")
                amine_reactant = any(
                    mol.HasSubstructMatch(amine_pattern) for mol in reactants if mol
                )

                # Check for amide in product
                amide_pattern = Chem.MolFromSmarts("[#6]-[C](=[O])[N]([#6])[#6]")
                product_has_amide = product.HasSubstructMatch(amide_pattern) if product else False

                # Check if this reaction is an amide coupling
                if (acid_reactant or ester_reactant) and amine_reactant and product_has_amide:
                    print(f"Found amide coupling at depth {depth}")
                    amide_couplings.append(depth)

        # Traverse children
        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    # Start traversal
    dfs_traverse(route)

    # Strategy is present if we found at least two amide couplings
    strategy_present = len(amide_couplings) >= 2

    print(f"Amide coupling fragment strategy detected: {strategy_present}")
    if strategy_present:
        print(f"Amide couplings occurred at depths: {amide_couplings}")

    return strategy_present
