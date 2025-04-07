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
    Detects synthesis routes that include an ester hydrolysis step.
    """
    found_ester_hydrolysis = False

    def dfs_traverse(node):
        nonlocal found_ester_hydrolysis

        if node["type"] == "reaction" and "metadata" in node and "rsmi" in node["metadata"]:
            rsmi = node["metadata"]["rsmi"]
            reactants = rsmi.split(">")[0].split(".")
            product = rsmi.split(">")[-1]

            reactant_mols = [Chem.MolFromSmiles(r) for r in reactants if r]
            product_mol = Chem.MolFromSmiles(product) if product else None

            if all(reactant_mols) and product_mol:
                # Ester pattern
                ester_pattern = Chem.MolFromSmarts("[#6]-[#8]-[#6](=[#8])")

                # Carboxylic acid pattern
                acid_pattern = Chem.MolFromSmarts("[#8;H1]-[#6](=[#8])")

                # Check for ester hydrolysis
                reactants_have_ester = any(
                    mol.HasSubstructMatch(ester_pattern) for mol in reactant_mols
                )
                product_has_acid = product_mol.HasSubstructMatch(acid_pattern)

                if reactants_have_ester and product_has_acid:
                    found_ester_hydrolysis = True
                    print("Found ester hydrolysis reaction")

        # Traverse children
        for child in node.get("children", []):
            dfs_traverse(child)

    # Start traversal
    dfs_traverse(route)

    print(f"Ester hydrolysis in sequence detected: {found_ester_hydrolysis}")
    return found_ester_hydrolysis
