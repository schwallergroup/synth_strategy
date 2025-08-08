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
    Detects if the synthesis route uses SNAr reaction for C-O bond formation.
    Looks for replacement of F with OR group on aromatic ring.
    """
    found = False

    def dfs_traverse(node):
        nonlocal found

        if node["type"] == "reaction" and "metadata" in node and "rsmi" in node["metadata"]:
            rsmi = node["metadata"]["rsmi"]
            reactants = rsmi.split(">")[0].split(".")
            product = rsmi.split(">")[-1]

            # Fluoro-aromatic pattern
            fluoro_pattern = Chem.MolFromSmarts("[c][F]")
            # Alkoxy pattern
            alkoxy_pattern = Chem.MolFromSmarts("[c][O][C]")

            reactant_mols = [Chem.MolFromSmiles(r) for r in reactants]
            product_mol = Chem.MolFromSmiles(product)

            has_fluoro = any(
                mol is not None and mol.HasSubstructMatch(fluoro_pattern) for mol in reactant_mols
            )
            has_alkoxy_in_product = product_mol is not None and product_mol.HasSubstructMatch(
                alkoxy_pattern
            )

            if has_fluoro and has_alkoxy_in_product:
                # Check if the alkoxy group replaced a fluorine
                # This is a simplification - a more robust implementation would check atom indices
                print("Found potential SNAr reaction")
                found = True

        for child in node.get("children", []):
            dfs_traverse(child)

    dfs_traverse(route)
    return found
