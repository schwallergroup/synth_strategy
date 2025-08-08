#!/bin/python

"""LM-defined function for strategy description."""

from rdkit.Chem import AllChem, rdFMCS
import copy
from collections import deque
import rdkit.Chem as Chem
from rdkit.Chem import rdMolDescriptors
from rdkit.Chem import rdChemReactions
from rdkit.Chem import AllChem
from rdkit.Chem import rdFMCS
import rdkit.Chem.rdFMCS
from rdkit.Chem import AllChem, Descriptors, rdMolDescriptors
from rdkit import Chem
from rdkit.Chem import Descriptors
from rdkit.Chem import AllChem, rdMolDescriptors
from rdkit.Chem import AllChem, Descriptors, Lipinski
from rdkit.Chem import rdmolops
import re
from rdkit.Chem.Scaffolds import MurckoScaffold
from rdkit.Chem import AllChem, Descriptors
import traceback
import rdkit
from collections import Counter


def main(route):
    """
    This function detects if the synthetic route involves reduction of carboxylic acid to alcohol.
    """
    carboxylic_acid_pattern = Chem.MolFromSmarts("C(=O)O")
    primary_alcohol_pattern = Chem.MolFromSmarts("CO")
    reduction_detected = False

    def dfs_traverse(node):
        nonlocal reduction_detected

        if node["type"] == "reaction":
            # Extract reactants and product
            rsmi = node["metadata"]["rsmi"]
            reactants_smiles = rsmi.split(">")[0].split(".")
            product_smiles = rsmi.split(">")[-1]

            # Check if any reactant contains carboxylic acid
            reactant_has_carboxylic_acid = False
            for r_smiles in reactants_smiles:
                try:
                    r_mol = Chem.MolFromSmiles(r_smiles)
                    if r_mol and r_mol.HasSubstructMatch(carboxylic_acid_pattern):
                        reactant_has_carboxylic_acid = True
                        break
                except:
                    continue

            # Check if product contains primary alcohol
            try:
                p_mol = Chem.MolFromSmiles(product_smiles)
                product_has_primary_alcohol = p_mol and p_mol.HasSubstructMatch(
                    primary_alcohol_pattern
                )
            except:
                product_has_primary_alcohol = False

            # If reactant has carboxylic acid and product has primary alcohol, it's a reduction
            if reactant_has_carboxylic_acid and product_has_primary_alcohol:
                print(f"Carboxylic acid reduction detected in reaction: {rsmi}")
                reduction_detected = True

        # Traverse children
        for child in node.get("children", []):
            dfs_traverse(child)

    # Start traversal
    dfs_traverse(route)
    return reduction_detected
