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
    This function detects if the synthetic route involves late-stage introduction
    of an amine group, particularly N,N-dimethylamine.
    """
    late_amine_introduction = False

    def dfs_traverse(node, depth=0):
        nonlocal late_amine_introduction

        if node["type"] == "reaction" and depth <= 1:  # Late stage (depth 0 or 1)
            rsmi = node["metadata"]["rsmi"]
            reactants_smiles = rsmi.split(">")[0]
            product_smiles = rsmi.split(">")[-1]

            try:
                # Look for dimethylamine or similar amine in reactants
                amine_pattern = Chem.MolFromSmarts("[#7]([#6])[#6]")

                for reactant in reactants_smiles.split("."):
                    reactant_mol = Chem.MolFromSmiles(reactant)
                    if reactant_mol and reactant_mol.HasSubstructMatch(amine_pattern):
                        # Check if amine is incorporated into product
                        product_mol = Chem.MolFromSmiles(product_smiles)
                        if product_mol and product_mol.HasSubstructMatch(amine_pattern):
                            print(f"Found late-stage amine introduction at depth {depth}")
                            late_amine_introduction = True
            except:
                pass

        # Traverse children
        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    dfs_traverse(route)
    return late_amine_introduction
