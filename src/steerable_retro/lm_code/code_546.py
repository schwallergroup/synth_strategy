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
    Detects if the synthesis route includes a nucleophilic aromatic substitution (SNAr)
    reaction where an amine replaces a halogen.
    """
    has_snar = False

    def dfs_traverse(node, depth=0):
        nonlocal has_snar

        if node["type"] == "reaction":
            rsmi = node["metadata"]["rsmi"]
            reactants_smiles = rsmi.split(">")[0].split(".")
            product_smiles = rsmi.split(">")[-1]

            # Look for chloro-aromatic and amine reactants
            chloro_aromatic_pattern = Chem.MolFromSmarts("c[Cl]")
            amine_pattern = Chem.MolFromSmarts("[NH2][C]")

            has_chloro_aromatic = False
            has_amine = False

            for reactant in reactants_smiles:
                mol = Chem.MolFromSmiles(reactant)
                if mol:
                    if mol.HasSubstructMatch(chloro_aromatic_pattern):
                        has_chloro_aromatic = True
                    if mol.HasSubstructMatch(amine_pattern):
                        has_amine = True

            # Check if product has C-N bond where chlorine was
            if has_chloro_aromatic and has_amine:
                product_mol = Chem.MolFromSmiles(product_smiles)
                if product_mol:
                    c_n_bond_pattern = Chem.MolFromSmarts("c[NH][C]")
                    if product_mol.HasSubstructMatch(c_n_bond_pattern):
                        has_snar = True
                        print(f"Detected SNAr amine introduction at depth {depth}")

        # Continue traversing
        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    # Start traversal
    dfs_traverse(route)
    return has_snar
