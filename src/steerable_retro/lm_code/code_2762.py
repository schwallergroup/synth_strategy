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
    This function detects a synthetic strategy involving the conversion of a
    Weinreb amide to a ketone.
    """
    # Initialize tracking variables
    has_weinreb_to_ketone = False

    def dfs_traverse(node):
        nonlocal has_weinreb_to_ketone

        if node["type"] == "reaction":
            # Extract reactants and product
            rsmi = node.get("metadata", {}).get("rsmi", "")
            if not rsmi:
                return

            reactants_smiles = rsmi.split(">")[0].split(".")
            product_smiles = rsmi.split(">")[-1]

            # Check for Weinreb amide in reactants
            weinreb_amide_pattern = Chem.MolFromSmarts("[#6]-[#8]-[#7](-[#6])-[#6](=[#8])-[#6]")
            # Check for ketone in product
            ketone_pattern = Chem.MolFromSmarts("[#6]-[#6](=[#8])-[#6]")

            has_weinreb = False
            for reactant in reactants_smiles:
                try:
                    mol = Chem.MolFromSmiles(reactant)
                    if mol and mol.HasSubstructMatch(weinreb_amide_pattern):
                        has_weinreb = True
                        break
                except:
                    continue

            has_ketone = False
            try:
                product_mol = Chem.MolFromSmiles(product_smiles)
                if product_mol and product_mol.HasSubstructMatch(ketone_pattern):
                    has_ketone = True
            except:
                pass

            if has_weinreb and has_ketone:
                print("Found Weinreb amide to ketone transformation")
                has_weinreb_to_ketone = True

        # Traverse children
        for child in node.get("children", []):
            dfs_traverse(child)

    # Start traversal
    dfs_traverse(route)

    return has_weinreb_to_ketone
