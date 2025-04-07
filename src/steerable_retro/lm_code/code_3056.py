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
    This function detects a synthetic strategy involving benzyl protection of hydroxyl groups,
    looking for at least one benzylation reaction.
    """
    has_benzylation = False

    def dfs_traverse(node):
        nonlocal has_benzylation

        if node["type"] == "reaction":
            # Extract reactants and product from reaction SMILES
            rsmi = node["metadata"]["rsmi"]
            reactants = rsmi.split(">")[0].split(".")
            product = rsmi.split(">")[-1]

            # Convert to RDKit molecules
            try:
                product_mol = Chem.MolFromSmiles(product)
                reactant_mols = [Chem.MolFromSmiles(r) for r in reactants if r.strip()]

                # Pattern for phenol in reactants
                phenol_pattern = Chem.MolFromSmarts("[OH]-[c]")

                # Pattern for benzyl ether in product
                benzyl_ether_pattern = Chem.MolFromSmarts("[O]-[CH2]-[c]")

                # Check for benzyl bromide or similar benzylating agent
                benzyl_halide_pattern = Chem.MolFromSmarts("c1ccccc1-[CH2]-[Br,Cl,I,O]")

                for reactant in reactant_mols:
                    if reactant and product_mol:
                        # Check if one reactant has phenol and another has benzyl halide
                        if reactant.HasSubstructMatch(phenol_pattern):
                            if any(
                                r and r.HasSubstructMatch(benzyl_halide_pattern)
                                for r in reactant_mols
                            ):
                                if product_mol.HasSubstructMatch(benzyl_ether_pattern):
                                    print(
                                        f"Benzyl protection detected at depth {node.get('depth', 'unknown')}"
                                    )
                                    has_benzylation = True
                                    break
            except:
                print("Error processing reaction SMILES")

        # Process children
        for child in node.get("children", []):
            dfs_traverse(child)

    # Start traversal
    dfs_traverse(route)

    return has_benzylation
