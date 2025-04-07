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
    Detects if the synthesis route involves multiple amide bond formations.
    """
    amide_formation_count = 0

    def dfs_traverse(node):
        nonlocal amide_formation_count

        if node["type"] == "reaction":
            # Extract reactants and product
            rsmi = node["metadata"]["rsmi"]
            reactants_smiles = rsmi.split(">")[0].split(".")
            product_smiles = rsmi.split(">")[-1]

            # Check if this is an amide formation
            reactant_mols = [Chem.MolFromSmiles(r) for r in reactants_smiles]
            product_mol = Chem.MolFromSmiles(product_smiles)

            # Look for amine in reactants
            amine_pattern = Chem.MolFromSmarts("[#7;H1,H2]")
            has_amine = any(
                mol.HasSubstructMatch(amine_pattern)
                for mol in reactant_mols
                if mol is not None
            )

            # Look for carbonyl-containing reactant (acid, acid chloride, etc.)
            carbonyl_pattern = Chem.MolFromSmarts("[#6](=[#8])[#8,#17]")
            has_carbonyl = any(
                mol.HasSubstructMatch(carbonyl_pattern)
                for mol in reactant_mols
                if mol is not None
            )

            # Look for amide in product
            amide_pattern = Chem.MolFromSmarts("[#6](=[#8])[#7]")
            has_amide_product = (
                product_mol is not None and product_mol.HasSubstructMatch(amide_pattern)
            )

            if has_amine and has_carbonyl and has_amide_product:
                amide_formation_count += 1
                print(f"Detected amide formation at reaction node")

        # Traverse children
        for child in node.get("children", []):
            dfs_traverse(child)

    # Start traversal from root
    dfs_traverse(route)
    print(f"Total amide formations detected: {amide_formation_count}")
    return amide_formation_count >= 2
