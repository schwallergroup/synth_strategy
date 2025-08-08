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
    This function detects an N-alkylation strategy, specifically
    alkylation of an indole or similar heterocycle nitrogen.
    """
    # Track if we found N-alkylation
    found_n_alkylation = False

    def dfs_traverse(node):
        nonlocal found_n_alkylation

        if node["type"] == "reaction":
            # Get reactants and product
            rsmi = node.get("metadata", {}).get("rsmi", "")
            if rsmi:
                reactants_smiles = rsmi.split(">")[0].split(".")
                product_smiles = rsmi.split(">")[-1]

                # Check for N-alkylation pattern
                nh_pattern = Chem.MolFromSmarts("[nH]")
                n_alkyl_pattern = Chem.MolFromSmarts("[n][C]")
                alkyl_halide_pattern = Chem.MolFromSmarts("[C][Br,Cl,I]")

                # Convert SMILES to molecules
                reactant_mols = [Chem.MolFromSmiles(r) for r in reactants_smiles if r]
                product_mol = Chem.MolFromSmiles(product_smiles) if product_smiles else None

                if product_mol and len(reactant_mols) >= 2:
                    # Check if any reactant has NH and any has alkyl halide
                    has_nh = any(mol and mol.HasSubstructMatch(nh_pattern) for mol in reactant_mols)
                    has_alkyl_halide = any(
                        mol and mol.HasSubstructMatch(alkyl_halide_pattern) for mol in reactant_mols
                    )

                    # Check if product has N-alkyl
                    has_n_alkyl = product_mol.HasSubstructMatch(n_alkyl_pattern)

                    if has_nh and has_alkyl_halide and has_n_alkyl:
                        found_n_alkylation = True
                        print("Found N-alkylation reaction")

        # Traverse children
        for child in node.get("children", []):
            dfs_traverse(child)

    # Start traversal
    dfs_traverse(route)

    return found_n_alkylation
