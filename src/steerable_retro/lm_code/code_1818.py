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
    Detects a strategy involving the formation of a quaternary carbon center.
    """
    found_quaternary_carbon_formation = False

    def dfs_traverse(node, depth=0):
        nonlocal found_quaternary_carbon_formation

        if node["type"] == "reaction":
            # Extract reactants and product
            rsmi = node["metadata"]["rsmi"]
            reactants_smiles = rsmi.split(">")[0].split(".")
            product_smiles = rsmi.split(">")[-1]

            # Convert to RDKit molecules
            reactant_mols = [Chem.MolFromSmiles(r) for r in reactants_smiles if r]
            product_mol = Chem.MolFromSmiles(product_smiles)

            if product_mol and all(r for r in reactant_mols):
                # Check for quaternary carbon in product
                quat_carbon_patt = Chem.MolFromSmarts("[#6D4](-[#6])(-[#6])(-[#6])(-[#6])")

                if product_mol.HasSubstructMatch(quat_carbon_patt):
                    # Check if quaternary carbon was formed in this step
                    reactants_have_quat_carbon = False
                    for r in reactant_mols:
                        if r.HasSubstructMatch(quat_carbon_patt):
                            reactants_have_quat_carbon = True
                            break

                    if not reactants_have_quat_carbon:
                        found_quaternary_carbon_formation = True
                        print(f"Found quaternary carbon formation at depth {depth}")

        # Traverse children
        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    # Start traversal
    dfs_traverse(route)

    return found_quaternary_carbon_formation
