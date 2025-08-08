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
    Detects if the synthesis includes the reduction of a cyclic amide (lactam) to a cyclic amine.
    """
    lactam_reduction_found = False

    def dfs_traverse(node, depth=0):
        nonlocal lactam_reduction_found

        if node["type"] == "reaction":
            rsmi = node["metadata"]["rsmi"]
            reactants_smiles = rsmi.split(">")[0].split(".")
            product_smiles = rsmi.split(">")[-1]

            # Check for lactam in reactants
            lactam_pattern = Chem.MolFromSmarts("[C]1[C][C][N][C](=[O])[C]1")

            reactants_have_lactam = False
            for reactant in reactants_smiles:
                try:
                    mol = Chem.MolFromSmiles(reactant)
                    if mol and mol.HasSubstructMatch(lactam_pattern):
                        reactants_have_lactam = True
                        break
                except:
                    continue

            # Check for cyclic amine in product
            cyclic_amine_pattern = Chem.MolFromSmarts("[C]1[C][C][N][C][C]1")

            product_has_cyclic_amine = False
            try:
                product_mol = Chem.MolFromSmiles(product_smiles)
                if product_mol and product_mol.HasSubstructMatch(cyclic_amine_pattern):
                    product_has_cyclic_amine = True
            except:
                pass

            if reactants_have_lactam and product_has_cyclic_amine:
                lactam_reduction_found = True
                print(f"Lactam reduction detected at depth {depth}")

        # Traverse children
        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    dfs_traverse(route)
    return lactam_reduction_found
