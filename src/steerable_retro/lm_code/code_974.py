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
    Detects a strategy involving the transformation of an ester to an aldehyde.
    """
    ester_to_aldehyde_found = False

    def dfs_traverse(node):
        nonlocal ester_to_aldehyde_found

        if (
            node.get("type") == "reaction"
            and "metadata" in node
            and "rsmi" in node.get("metadata", {})
        ):
            rsmi = node["metadata"]["rsmi"]
            reactants_part = rsmi.split(">")[0]
            product_part = rsmi.split(">")[-1]

            # Check for ester in reactants
            ester_pattern = Chem.MolFromSmarts("[#6]-[#8]-[#6](=[#8])")

            # Check for aldehyde in product
            aldehyde_pattern = Chem.MolFromSmarts("[#6;H1]=[#8]")

            reactants = [Chem.MolFromSmiles(r) for r in reactants_part.split(".") if r]
            product = Chem.MolFromSmiles(product_part) if product_part else None

            if product and all(r for r in reactants):
                has_ester = any(r.HasSubstructMatch(ester_pattern) for r in reactants if r)
                has_aldehyde = product.HasSubstructMatch(aldehyde_pattern)

                if has_ester and has_aldehyde:
                    print("Ester to aldehyde transformation detected")
                    ester_to_aldehyde_found = True

        for child in node.get("children", []):
            dfs_traverse(child)

    dfs_traverse(route)

    print(f"Ester to aldehyde transformation found: {ester_to_aldehyde_found}")

    return ester_to_aldehyde_found
