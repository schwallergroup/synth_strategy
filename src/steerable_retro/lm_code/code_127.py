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
    This function detects a synthetic strategy involving the transformation
    of a nitrile functional group to a carboxylic acid.
    """
    # Track if we found nitrile to acid transformation
    found_nitrile_to_acid = False

    def dfs_traverse(node, depth=0):
        nonlocal found_nitrile_to_acid

        if node["type"] == "reaction":
            if "metadata" in node and "rsmi" in node["metadata"]:
                rsmi = node["metadata"]["rsmi"]
                reactants_smiles = rsmi.split(">")[0]
                product_smiles = rsmi.split(">")[-1]

                reactants = [Chem.MolFromSmiles(r) for r in reactants_smiles.split(".") if r]
                product = Chem.MolFromSmiles(product_smiles) if product_smiles else None

                if product and reactants:
                    # Check for nitrile in reactants
                    nitrile_pattern = Chem.MolFromSmarts("[CX2]#[NX1]")
                    reactant_has_nitrile = any(
                        r.HasSubstructMatch(nitrile_pattern) for r in reactants if r
                    )

                    # Check for carboxylic acid in product
                    acid_pattern = Chem.MolFromSmarts("[CX3](=O)[OX2H1]")
                    product_has_acid = product.HasSubstructMatch(acid_pattern)

                    if reactant_has_nitrile and product_has_acid:
                        print(f"Found nitrile to carboxylic acid transformation at depth {depth}")
                        found_nitrile_to_acid = True

        # Traverse children
        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    # Start traversal
    dfs_traverse(route)

    print(f"Nitrile to acid transformation detected: {found_nitrile_to_acid}")
    return found_nitrile_to_acid
