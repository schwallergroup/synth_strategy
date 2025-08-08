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
    Detects benzylic bromination as a key functional group transformation.
    """
    benzylic_bromination_found = False

    def dfs_traverse(node):
        nonlocal benzylic_bromination_found

        if node["type"] == "reaction":
            rsmi = node["metadata"]["rsmi"]
            reactants_smiles = rsmi.split(">")[0].split(".")
            product_smiles = rsmi.split(">")[-1]

            # Check for benzylic methyl in reactants
            benzylic_methyl_in_reactants = False
            for r in reactants_smiles:
                r_mol = Chem.MolFromSmiles(r)
                if r_mol:
                    benzylic_methyl_pattern = Chem.MolFromSmarts("c[CH3]")
                    if r_mol.HasSubstructMatch(benzylic_methyl_pattern):
                        benzylic_methyl_in_reactants = True
                        break

            # Check for benzylic bromide in product
            if benzylic_methyl_in_reactants:
                product_mol = Chem.MolFromSmiles(product_smiles)
                if product_mol:
                    benzylic_bromide_pattern = Chem.MolFromSmarts("c[CH2][Br]")
                    if product_mol.HasSubstructMatch(benzylic_bromide_pattern):
                        benzylic_bromination_found = True
                        print(f"Found benzylic bromination: {rsmi}")

        # Traverse children
        for child in node.get("children", []):
            dfs_traverse(child)

    # Start traversal
    dfs_traverse(route)

    return benzylic_bromination_found
