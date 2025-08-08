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
    Detects if the synthesis uses chloroacetyl chloride as a linker/handle
    for subsequent substitution reactions.
    """
    has_chloroacetyl_chloride = False
    has_chloroacetamide_intermediate = False
    has_substitution_product = False

    def dfs_traverse(node):
        nonlocal has_chloroacetyl_chloride, has_chloroacetamide_intermediate, has_substitution_product

        if node["type"] == "reaction":
            rsmi = node.get("metadata", {}).get("rsmi", "")
            if not rsmi:
                return

            reactants_smiles = rsmi.split(">")[0].split(".")
            product_smiles = rsmi.split(">")[-1]

            # Check for chloroacetyl chloride
            chloroacetyl_pattern = Chem.MolFromSmarts("Cl[C](=[O])[CH2][Cl]")
            for reactant in reactants_smiles:
                try:
                    mol = Chem.MolFromSmiles(reactant)
                    if mol and mol.HasSubstructMatch(chloroacetyl_pattern):
                        has_chloroacetyl_chloride = True
                        print("Found chloroacetyl chloride")
                except:
                    continue

            # Check for chloroacetamide intermediate
            chloroacetamide_pattern = Chem.MolFromSmarts("[NH][C](=[O])[CH2][Cl]")
            try:
                product_mol = Chem.MolFromSmiles(product_smiles)
                if product_mol and product_mol.HasSubstructMatch(chloroacetamide_pattern):
                    has_chloroacetamide_intermediate = True
                    print("Found chloroacetamide intermediate")
            except:
                pass

            # Check for substitution product
            substitution_product_pattern = Chem.MolFromSmarts("[NH][C](=[O])[CH2][N]")
            try:
                product_mol = Chem.MolFromSmiles(product_smiles)
                if product_mol and product_mol.HasSubstructMatch(substitution_product_pattern):
                    has_substitution_product = True
                    print("Found substitution product")
            except:
                pass

        # Traverse children
        for child in node.get("children", []):
            dfs_traverse(child)

    # Start traversal from the root
    dfs_traverse(route)

    # Strategy is present if all three conditions are met
    return (
        has_chloroacetyl_chloride and has_chloroacetamide_intermediate and has_substitution_product
    )
