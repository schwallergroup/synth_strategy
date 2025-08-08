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
    Detects if the synthesis route involves a diaryl ether disconnection as a key step.
    """
    diaryl_ether_disconnection_found = False

    def dfs_traverse(node):
        nonlocal diaryl_ether_disconnection_found

        if node["type"] == "reaction":
            # Extract reactants and product
            rsmi = node["metadata"]["rsmi"]
            reactants_smiles = rsmi.split(">")[0].split(".")
            product_smiles = rsmi.split(">")[-1]

            # Check for diaryl ether disconnection
            reactant_mols = [Chem.MolFromSmiles(smi) for smi in reactants_smiles]
            product_mol = Chem.MolFromSmiles(product_smiles)

            if product_mol and all(r for r in reactant_mols if r):
                # Check for phenol pattern in reactants
                phenol_pattern = Chem.MolFromSmarts("[c][OH]")
                # Check for diaryl ether pattern in product
                diaryl_ether_pattern = Chem.MolFromSmarts("[c][O][c]")

                has_phenol = any(r.HasSubstructMatch(phenol_pattern) for r in reactant_mols if r)
                has_diaryl_ether = product_mol.HasSubstructMatch(diaryl_ether_pattern)

                if has_phenol and has_diaryl_ether:
                    diaryl_ether_disconnection_found = True
                    print("Diaryl ether disconnection found")

        # Traverse children
        for child in node.get("children", []):
            dfs_traverse(child)

    dfs_traverse(route)

    print(f"Diaryl ether disconnection strategy: {diaryl_ether_disconnection_found}")
    return diaryl_ether_disconnection_found
