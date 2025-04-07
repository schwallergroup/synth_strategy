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
    This function detects a late-stage O-alkylation, specifically
    the formation of a methoxymethyl ether in the final synthetic step.
    """
    o_alkylation_detected = False
    depth_of_alkylation = float("inf")

    def dfs_traverse(node, depth=0):
        nonlocal o_alkylation_detected, depth_of_alkylation

        if node["type"] == "reaction":
            if "rsmi" in node.get("metadata", {}):
                rsmi = node["metadata"]["rsmi"]
                reactants = rsmi.split(">")[0].split(".")
                product = rsmi.split(">")[-1]

                # Check for alcohol in reactants
                alcohol_pattern = Chem.MolFromSmarts("[OH]")

                # Check for chloromethyl ether in reactants
                chloromethyl_ether_pattern = Chem.MolFromSmarts("Cl[CH2]O[CH3]")

                # Check for methoxymethyl ether in product
                methoxymethyl_ether_pattern = Chem.MolFromSmarts("[#6]O[CH2]O[CH3]")

                # Convert SMILES to RDKit molecules
                reactant_mols = [Chem.MolFromSmiles(r) for r in reactants if r]
                product_mol = Chem.MolFromSmiles(product) if product else None

                if (
                    product_mol
                    and any(mol and mol.HasSubstructMatch(alcohol_pattern) for mol in reactant_mols)
                    and any(
                        mol and mol.HasSubstructMatch(chloromethyl_ether_pattern)
                        for mol in reactant_mols
                    )
                    and product_mol.HasSubstructMatch(methoxymethyl_ether_pattern)
                ):
                    print(f"Detected O-alkylation at depth {depth}")
                    o_alkylation_detected = True
                    depth_of_alkylation = min(depth_of_alkylation, depth)

        # Traverse children
        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    # Start traversal
    dfs_traverse(route)

    # Consider it late-stage if it occurs at depth 0 or 1
    return o_alkylation_detected and depth_of_alkylation <= 1
