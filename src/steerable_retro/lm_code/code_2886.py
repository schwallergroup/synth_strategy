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
    Detects a linear synthesis route with multiple amide bond formations.
    """
    amide_formation_count = 0
    is_linear = True

    def dfs_traverse(node):
        nonlocal amide_formation_count, is_linear

        if node["type"] == "reaction":
            # Extract reactants and product
            rsmi = node["metadata"]["rsmi"]
            reactants_smiles = rsmi.split(">")[0].split(".")
            product_smiles = rsmi.split(">")[-1].split(".")

            # Check if this is a linear step (exactly one product)
            if len(product_smiles) != 1:
                is_linear = False

            # Check for amide formation
            has_carboxylic_acid = False
            has_amine = False

            for reactant in reactants_smiles:
                reactant_mol = Chem.MolFromSmiles(reactant)
                if reactant_mol:
                    if reactant_mol.HasSubstructMatch(Chem.MolFromSmarts("[C](=[O])[OH]")):
                        has_carboxylic_acid = True
                    if reactant_mol.HasSubstructMatch(Chem.MolFromSmarts("[NH2]")):
                        has_amine = True

            product_mol = Chem.MolFromSmiles(product_smiles[0])
            has_amide_product = False
            if product_mol and product_mol.HasSubstructMatch(Chem.MolFromSmarts("[C](=[O])[NH]")):
                has_amide_product = True

            if has_carboxylic_acid and has_amine and has_amide_product:
                amide_formation_count += 1

        # Process children
        for child in node.get("children", []):
            dfs_traverse(child)

    # Start traversal
    dfs_traverse(route)

    if is_linear and amide_formation_count >= 2:
        print(f"Found linear synthesis with {amide_formation_count} amide formations")
        return True
    else:
        print(
            f"Not a linear synthesis with multiple amide formations: linear={is_linear}, amide_count={amide_formation_count}"
        )
        return False
