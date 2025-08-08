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
    Detects a strategy involving Boc deprotection of an amine.
    """
    boc_deprotection_found = False

    def dfs_traverse(node):
        nonlocal boc_deprotection_found

        if (
            node.get("type") == "reaction"
            and "metadata" in node
            and "rsmi" in node.get("metadata", {})
        ):
            rsmi = node["metadata"]["rsmi"]
            reactants_part = rsmi.split(">")[0]
            product_part = rsmi.split(">")[-1]

            # Check for Boc-protected amine in reactants
            boc_amine_pattern = Chem.MolFromSmarts("[#6]-[#8]-[#6](=[#8])-[#7]")

            # Check for free amine in product
            free_amine_pattern = Chem.MolFromSmarts("[#7;H2]-[#6]")

            reactants = [Chem.MolFromSmiles(r) for r in reactants_part.split(".") if r]
            product = Chem.MolFromSmiles(product_part) if product_part else None

            if product and all(r for r in reactants):
                has_boc_amine = any(r.HasSubstructMatch(boc_amine_pattern) for r in reactants if r)
                has_free_amine = product.HasSubstructMatch(free_amine_pattern)

                if has_boc_amine and has_free_amine:
                    print("Boc deprotection detected")
                    boc_deprotection_found = True

        for child in node.get("children", []):
            dfs_traverse(child)

    dfs_traverse(route)

    print(f"Boc deprotection found: {boc_deprotection_found}")

    return boc_deprotection_found
