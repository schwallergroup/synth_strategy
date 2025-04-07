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
    This function detects if the synthetic route involves nucleophilic aromatic substitution
    of a chloropyrimidine with an amine.
    """
    nas_detected = False

    def dfs_traverse(node):
        nonlocal nas_detected

        if node["type"] == "reaction":
            # Extract reactants and product
            rsmi = node["metadata"]["rsmi"]
            reactants_smiles = rsmi.split(">")[0].split(".")
            product_smiles = rsmi.split(">")[-1]

            # Convert to RDKit molecules
            reactants = [Chem.MolFromSmiles(r) for r in reactants_smiles]
            product = Chem.MolFromSmiles(product_smiles)

            # Check for chloropyrimidine in reactants
            chloropyrimidine_pattern = Chem.MolFromSmarts("[#17]-c1ncncc1")

            # Check for aminopyrimidine in product
            aminopyrimidine_pattern = Chem.MolFromSmarts(
                "[#7;!H0,!$(N-[#6]=O)]-c1ncncc1"
            )

            # Check for amine in reactants
            amine_pattern = Chem.MolFromSmarts("[#7;!$(N=*);!$(NC=O)]")

            reactant_has_chloropyrimidine = any(
                r is not None and r.HasSubstructMatch(chloropyrimidine_pattern)
                for r in reactants
            )
            reactant_has_amine = any(
                r is not None and r.HasSubstructMatch(amine_pattern) for r in reactants
            )
            product_has_aminopyrimidine = (
                product is not None
                and product.HasSubstructMatch(aminopyrimidine_pattern)
            )

            if (
                reactant_has_chloropyrimidine
                and reactant_has_amine
                and product_has_aminopyrimidine
            ):
                print(
                    "Detected nucleophilic aromatic substitution of chloropyrimidine with amine"
                )
                nas_detected = True

        # Traverse children
        for child in node.get("children", []):
            dfs_traverse(child)

    # Start traversal
    dfs_traverse(route)
    return nas_detected
