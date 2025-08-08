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
    This function detects a strategy involving sequential aromatic nucleophilic substitutions
    where amine nucleophiles displace halogens on aromatic rings.
    """
    aromatic_subst_count = 0

    def dfs_traverse(node):
        nonlocal aromatic_subst_count

        if node["type"] == "reaction":
            # Extract reactants and product
            rsmi = node["metadata"]["rsmi"]
            reactants = rsmi.split(">")[0].split(".")
            product = rsmi.split(">")[-1]

            # Check for aromatic nucleophilic substitution pattern
            # Look for halogen in reactants and amine-aromatic bond in product
            halogen_pattern = Chem.MolFromSmarts("[c][F,Cl,Br,I]")
            amine_pattern = Chem.MolFromSmarts("[c][NH][c,C]")

            halogen_found = False
            for reactant in reactants:
                mol = Chem.MolFromSmiles(reactant)
                if mol and mol.HasSubstructMatch(halogen_pattern):
                    halogen_found = True
                    break

            product_mol = Chem.MolFromSmiles(product)
            amine_bond_found = product_mol and product_mol.HasSubstructMatch(amine_pattern)

            if halogen_found and amine_bond_found:
                aromatic_subst_count += 1
                print(
                    f"Found aromatic nucleophilic substitution at depth {node.get('depth', 'unknown')}"
                )

        # Traverse children
        for child in node.get("children", []):
            dfs_traverse(child)

    # Start traversal
    dfs_traverse(route)

    # Return True if we found at least 2 sequential aromatic nucleophilic substitutions
    result = aromatic_subst_count >= 2
    print(f"Sequential aromatic nucleophilic substitution strategy detected: {result}")
    return result
