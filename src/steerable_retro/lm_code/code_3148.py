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
    Detects a strategy involving exchange of Boc protection group to ethoxycarbonyl.
    """
    found_protection_exchange = False

    def dfs_traverse(node):
        nonlocal found_protection_exchange

        if node["type"] == "reaction":
            # Extract reactants and product
            rsmi = node["metadata"]["rsmi"]
            reactants_smiles = rsmi.split(">")[0].split(".")
            product_smiles = rsmi.split(">")[-1]

            # Convert to RDKit molecules
            reactants = [Chem.MolFromSmiles(r) for r in reactants_smiles]
            product = Chem.MolFromSmiles(product_smiles)

            if product and all(r for r in reactants):
                # Check for Boc group in reactants
                boc_pattern = Chem.MolFromSmarts("CC(C)(C)OC(=O)")
                # Check for ethoxycarbonyl in product
                ethoxycarbonyl_pattern = Chem.MolFromSmarts("CCOC(=O)")

                reactants_with_boc = any(r.HasSubstructMatch(boc_pattern) for r in reactants if r)
                product_with_ethoxycarbonyl = product.HasSubstructMatch(ethoxycarbonyl_pattern)

                if reactants_with_boc and product_with_ethoxycarbonyl:
                    print("Found protection group exchange (Boc to ethoxycarbonyl)")
                    found_protection_exchange = True

        # Process children
        for child in node.get("children", []):
            dfs_traverse(child)

    # Start traversal
    dfs_traverse(route)

    return found_protection_exchange
