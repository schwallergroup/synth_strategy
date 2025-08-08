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
    Detects if the route involves oxidation of an amine (NH2) to a nitro group (NO2).
    """
    amine_to_nitro_found = False

    def dfs_traverse(node):
        nonlocal amine_to_nitro_found

        if node["type"] == "reaction":
            rsmi = node["metadata"]["rsmi"]
            reactants = rsmi.split(">")[0].split(".")
            product = rsmi.split(">")[-1]

            # Check for amine in reactants
            amine_pattern = Chem.MolFromSmarts("[NH2]")

            # Check for nitro in products
            nitro_pattern = Chem.MolFromSmarts("[N+](=[O])[O-]")

            reactants_have_amine = any(
                Chem.MolFromSmiles(r) and Chem.MolFromSmiles(r).HasSubstructMatch(amine_pattern)
                for r in reactants
                if Chem.MolFromSmiles(r)
            )

            product_has_nitro = Chem.MolFromSmiles(product) and Chem.MolFromSmiles(
                product
            ).HasSubstructMatch(nitro_pattern)

            if reactants_have_amine and product_has_nitro:
                print("Amine to nitro oxidation detected")
                amine_to_nitro_found = True

        # Process children
        for child in node.get("children", []):
            dfs_traverse(child)

    dfs_traverse(route)

    return amine_to_nitro_found
