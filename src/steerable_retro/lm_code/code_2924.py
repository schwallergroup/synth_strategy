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
    This function detects thiazole ring formation in the synthetic route.
    """
    thiazole_pattern = Chem.MolFromSmarts("c1nc([#6])sc1")
    thiazole_formed = False

    def dfs_traverse(node):
        nonlocal thiazole_formed

        if node["type"] == "reaction":
            if "rsmi" in node.get("metadata", {}):
                rsmi = node["metadata"]["rsmi"]
                reactants_smiles = rsmi.split(">")[0]
                product_smiles = rsmi.split(">")[-1]

                try:
                    reactants_mol = Chem.MolFromSmiles(reactants_smiles)
                    product_mol = Chem.MolFromSmiles(product_smiles)

                    if reactants_mol and product_mol:
                        # Check if thiazole is in product but not in reactants
                        if product_mol.HasSubstructMatch(
                            thiazole_pattern
                        ) and not reactants_mol.HasSubstructMatch(thiazole_pattern):
                            thiazole_formed = True
                            print(f"Thiazole formation detected in reaction: {rsmi}")
                except Exception as e:
                    print(f"Error processing reaction SMILES: {e}")

        for child in node.get("children", []):
            dfs_traverse(child)

    dfs_traverse(route)
    return thiazole_formed
