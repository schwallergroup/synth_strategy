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
    Detects installation of difluoromethoxy (OCF2H) group in the synthesis.
    """
    found_difluoromethoxy_installation = False

    def dfs_traverse(node):
        nonlocal found_difluoromethoxy_installation

        if node["type"] == "reaction" and node.get("metadata", {}).get("rsmi"):
            rsmi = node["metadata"]["rsmi"]
            reactants_smiles = rsmi.split(">")[0].split(".")
            product_smiles = rsmi.split(">")[-1]

            try:
                reactants = [Chem.MolFromSmiles(r) for r in reactants_smiles]
                product = Chem.MolFromSmiles(product_smiles)

                if product:
                    # Check for difluoromethoxy pattern in product
                    difluoromethoxy_pattern = Chem.MolFromSmarts("[#8][#6]([#9])[#9]")
                    if product.HasSubstructMatch(difluoromethoxy_pattern):
                        # Check if any reactant has a phenol or similar group
                        phenol_pattern = Chem.MolFromSmarts("c[OH]")
                        has_phenol = any(
                            mol.HasSubstructMatch(phenol_pattern) for mol in reactants if mol
                        )

                        if has_phenol:
                            print("Found difluoromethoxy installation")
                            found_difluoromethoxy_installation = True
            except:
                print("Error processing reaction SMILES for difluoromethoxy detection")

        # Process children
        for child in node.get("children", []):
            dfs_traverse(child)

    # Start traversal
    dfs_traverse(route)
    return found_difluoromethoxy_installation
