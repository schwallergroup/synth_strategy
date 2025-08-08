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
    This function detects if the synthetic route involves formation of an acid chloride
    from a carboxylic acid.
    """
    acid_pattern = Chem.MolFromSmarts("C(=O)O")
    acid_chloride_pattern = Chem.MolFromSmarts("C(=O)Cl")

    acid_to_acid_chloride = False

    def dfs_traverse(node):
        nonlocal acid_to_acid_chloride

        if node["type"] == "reaction":
            if "rsmi" in node.get("metadata", {}):
                rsmi = node["metadata"]["rsmi"]
                reactants_smiles = rsmi.split(">")[0].split(".")
                product_smiles = rsmi.split(">")[-1]

                try:
                    # Check for acid in reactants
                    reactants_have_acid = any(
                        Chem.MolFromSmiles(r)
                        and Chem.MolFromSmiles(r).HasSubstructMatch(acid_pattern)
                        for r in reactants_smiles
                        if Chem.MolFromSmiles(r)
                    )

                    # Check for acid chloride in product
                    product_mol = Chem.MolFromSmiles(product_smiles)
                    if (
                        reactants_have_acid
                        and product_mol
                        and product_mol.HasSubstructMatch(acid_chloride_pattern)
                    ):
                        acid_to_acid_chloride = True
                        print("Detected carboxylic acid to acid chloride transformation")
                except:
                    pass

        for child in node.get("children", []):
            dfs_traverse(child)

    dfs_traverse(route)
    return acid_to_acid_chloride
