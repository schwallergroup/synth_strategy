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
    This function detects if the synthesis involves reduction of a nitrile to a primary amine.
    """
    result = False

    def dfs_traverse(node):
        nonlocal result

        if node["type"] == "reaction" and "metadata" in node and "rsmi" in node["metadata"]:
            rsmi = node["metadata"]["rsmi"]
            reactants_part = rsmi.split(">")[0]
            product_part = rsmi.split(">")[-1]

            # Check for nitrile pattern in reactants and amine in products
            nitrile_pattern = Chem.MolFromSmarts("[C]#[N]")
            amine_pattern = Chem.MolFromSmarts("[CH2][NH2]")

            reactant_mol = Chem.MolFromSmiles(reactants_part)
            product_mol = Chem.MolFromSmiles(product_part)

            if reactant_mol and product_mol:
                has_nitrile = reactant_mol.HasSubstructMatch(nitrile_pattern)
                has_amine = product_mol.HasSubstructMatch(amine_pattern)

                # Check if nitrile is reduced to amine
                if has_nitrile and has_amine and not product_mol.HasSubstructMatch(nitrile_pattern):
                    print("Found nitrile reduction to amine")
                    result = True

        # Continue traversing
        for child in node.get("children", []):
            dfs_traverse(child)

    # Start traversal from the root
    dfs_traverse(route)
    return result
