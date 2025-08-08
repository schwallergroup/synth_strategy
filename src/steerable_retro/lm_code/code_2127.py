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
    This function detects nucleophilic aromatic substitution reactions
    involving fluorinated aromatics to form diaryl ethers.
    """
    snar_found = False

    def dfs_traverse(node):
        nonlocal snar_found

        if node["type"] == "reaction":
            rsmi = node["metadata"]["rsmi"]
            reactants_smiles = rsmi.split(">")[0].split(".")
            product_smiles = rsmi.split(">")[-1]

            reactant_mols = [Chem.MolFromSmiles(r) for r in reactants_smiles]
            product_mol = Chem.MolFromSmiles(product_smiles)

            # Check for fluorinated aromatic
            fluoro_aromatic_pattern = Chem.MolFromSmarts("cF")
            # Check for phenol
            phenol_pattern = Chem.MolFromSmarts("c[OH]")
            # Check for diaryl ether
            diaryl_ether_pattern = Chem.MolFromSmarts("cOc")

            has_fluoro_aromatic = any(
                mol and mol.HasSubstructMatch(fluoro_aromatic_pattern) for mol in reactant_mols
            )
            has_phenol = any(mol and mol.HasSubstructMatch(phenol_pattern) for mol in reactant_mols)
            has_diaryl_ether = product_mol and product_mol.HasSubstructMatch(diaryl_ether_pattern)

            if has_fluoro_aromatic and has_phenol and has_diaryl_ether:
                snar_found = True
                print(f"Detected nucleophilic aromatic substitution: {rsmi}")

        # Traverse children
        for child in node.get("children", []):
            dfs_traverse(child)

    # Start traversal
    dfs_traverse(route)

    return snar_found
