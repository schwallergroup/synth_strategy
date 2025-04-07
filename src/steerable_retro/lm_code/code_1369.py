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
    This function detects if the synthetic route involves a hydroxyl → chloro → amino
    transformation sequence on a pyrimidine ring.
    """
    # Track if we've seen each transformation
    hydroxyl_to_chloro = False
    chloro_to_amino = False

    def dfs_traverse(node):
        nonlocal hydroxyl_to_chloro, chloro_to_amino

        if node["type"] == "reaction":
            # Extract reactants and product
            rsmi = node["metadata"]["rsmi"]
            reactants_smiles = rsmi.split(">")[0].split(".")
            product_smiles = rsmi.split(">")[-1]

            # Convert to RDKit molecules
            reactants = [Chem.MolFromSmiles(r) for r in reactants_smiles]
            product = Chem.MolFromSmiles(product_smiles)

            # Check for hydroxyl to chloro transformation on pyrimidine
            hydroxypyrimidine = Chem.MolFromSmarts("[#8H1]-c1ncncc1")
            chloropyrimidine = Chem.MolFromSmarts("[#17]-c1ncncc1")

            reactant_has_hydroxypyrimidine = any(
                r is not None and r.HasSubstructMatch(hydroxypyrimidine) for r in reactants
            )
            product_has_chloropyrimidine = product is not None and product.HasSubstructMatch(
                chloropyrimidine
            )

            if reactant_has_hydroxypyrimidine and product_has_chloropyrimidine:
                print("Detected hydroxyl to chloro transformation on pyrimidine")
                hydroxyl_to_chloro = True

            # Check for chloro to amino transformation on pyrimidine
            aminopyrimidine = Chem.MolFromSmarts("[#7;!H0,!$(N-[#6]=O)]-c1ncncc1")

            reactant_has_chloropyrimidine = any(
                r is not None and r.HasSubstructMatch(chloropyrimidine) for r in reactants
            )
            product_has_aminopyrimidine = product is not None and product.HasSubstructMatch(
                aminopyrimidine
            )

            if reactant_has_chloropyrimidine and product_has_aminopyrimidine:
                print("Detected chloro to amino transformation on pyrimidine")
                chloro_to_amino = True

        # Traverse children
        for child in node.get("children", []):
            dfs_traverse(child)

    # Start traversal
    dfs_traverse(route)
    return hydroxyl_to_chloro and chloro_to_amino
