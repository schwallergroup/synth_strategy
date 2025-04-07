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
    This function detects if the synthetic route involves pyrimidine ring formation
    from acyclic precursors.
    """
    pyrimidine_formed = False

    def dfs_traverse(node):
        nonlocal pyrimidine_formed

        if node["type"] == "reaction":
            # Extract reactants and product
            rsmi = node["metadata"]["rsmi"]
            reactants_smiles = rsmi.split(">")[0].split(".")
            product_smiles = rsmi.split(">")[-1]

            # Convert to RDKit molecules
            reactants = [Chem.MolFromSmiles(r) for r in reactants_smiles]
            product = Chem.MolFromSmiles(product_smiles)

            # Check if any reactant has a pyrimidine ring
            reactant_has_pyrimidine = False
            for r in reactants:
                if r is not None and r.HasSubstructMatch(Chem.MolFromSmarts("c1ncncc1")):
                    reactant_has_pyrimidine = True
                    break

            # Check if product has a pyrimidine ring
            product_has_pyrimidine = False
            if product is not None and product.HasSubstructMatch(Chem.MolFromSmarts("c1ncncc1")):
                product_has_pyrimidine = True

            # If product has pyrimidine but reactants don't, it's a pyrimidine formation
            if product_has_pyrimidine and not reactant_has_pyrimidine:
                print("Detected pyrimidine ring formation")
                pyrimidine_formed = True

        # Traverse children
        for child in node.get("children", []):
            dfs_traverse(child)

    # Start traversal
    dfs_traverse(route)
    return pyrimidine_formed
