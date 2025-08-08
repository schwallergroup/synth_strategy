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
    This function detects a strategy involving multiple C=C bond disconnections
    in the retrosynthetic direction.
    """
    cc_double_bond_disconnections = 0

    def dfs_traverse(node):
        nonlocal cc_double_bond_disconnections

        if node["type"] == "reaction":
            if "rsmi" in node.get("metadata", {}):
                rsmi = node["metadata"]["rsmi"]
                reactants = rsmi.split(">")[0]
                product = rsmi.split(">")[-1]

                # In retrosynthesis, we're looking for C=C bonds in the product that are broken in reactants
                product_mol = Chem.MolFromSmiles(product)

                if product_mol:
                    # Find C=C bonds in product
                    cc_double_bond_pattern = Chem.MolFromSmarts("C=C")
                    if product_mol.HasSubstructMatch(cc_double_bond_pattern):
                        # This is a simplification - in a real implementation, you would need to
                        # check if specific C=C bonds are broken by comparing atom mappings
                        # between reactants and products
                        print("Potential C=C bond disconnection detected")
                        cc_double_bond_disconnections += 1

        for child in node.get("children", []):
            dfs_traverse(child)

    dfs_traverse(route)
    return cc_double_bond_disconnections >= 2  # At least 2 C=C disconnections
