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
    This function detects a synthetic strategy that involves multiple amide bond
    formations or disconnections.
    """
    amide_reactions_count = 0

    def dfs_traverse(node):
        nonlocal amide_reactions_count

        if node["type"] == "reaction":
            # Get reactants and product
            rsmi = node["metadata"]["rsmi"]
            reactants_smiles = rsmi.split(">")[0].split(".")
            product_smiles = rsmi.split(">")[-1]

            reactant_mols = [Chem.MolFromSmiles(r) for r in reactants_smiles]
            product_mol = Chem.MolFromSmiles(product_smiles)

            # Amide pattern
            amide_pattern = Chem.MolFromSmarts("[#6]-[#7]-[#6](=[#8])-[#6]")

            # Check for amide formation: product has amide but reactants don't
            if product_mol and product_mol.HasSubstructMatch(amide_pattern):
                reactants_with_amide = any(
                    mol and mol.HasSubstructMatch(amide_pattern) for mol in reactant_mols
                )
                if not reactants_with_amide:
                    print("Found amide formation")
                    amide_reactions_count += 1

            # Check for amide disconnection: reactants have amide but product doesn't
            reactants_with_amide = any(
                mol and mol.HasSubstructMatch(amide_pattern) for mol in reactant_mols
            )
            if reactants_with_amide and (
                not product_mol or not product_mol.HasSubstructMatch(amide_pattern)
            ):
                print("Found amide disconnection")
                amide_reactions_count += 1

        # Traverse children
        for child in node.get("children", []):
            dfs_traverse(child)

    # Start traversal
    dfs_traverse(route)

    # Return True if at least 2 amide reactions are found
    return amide_reactions_count >= 2
