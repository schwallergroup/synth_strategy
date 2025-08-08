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
    Detects if the synthesis uses a chloroacetamide intermediate as a reactive handle
    for late-stage incorporation of a piperazine fragment through nucleophilic substitution.
    """
    # Initialize tracking variables
    has_chloroacetamide_intermediate = False
    has_piperazine_nucleophile = False
    has_final_piperazine_product = False

    def dfs_traverse(node):
        nonlocal has_chloroacetamide_intermediate, has_piperazine_nucleophile, has_final_piperazine_product

        if node["type"] == "reaction":
            # Extract reactants and product from reaction SMILES
            rsmi = node.get("metadata", {}).get("rsmi", "")
            if not rsmi:
                return

            reactants_smiles = rsmi.split(">")[0].split(".")
            product_smiles = rsmi.split(">")[-1]

            # Check for chloroacetamide intermediate
            chloroacetamide_pattern = Chem.MolFromSmarts("[NH][C](=[O])[CH2][Cl]")
            for reactant in reactants_smiles:
                try:
                    mol = Chem.MolFromSmiles(reactant)
                    if mol and mol.HasSubstructMatch(chloroacetamide_pattern):
                        has_chloroacetamide_intermediate = True
                        print("Found chloroacetamide intermediate")
                except:
                    continue

            # Check for piperazine nucleophile
            piperazine_pattern = Chem.MolFromSmarts("[CH3][N]1[CH2][CH2][N][CH2][CH2]1")
            for reactant in reactants_smiles:
                try:
                    mol = Chem.MolFromSmiles(reactant)
                    if mol and mol.HasSubstructMatch(piperazine_pattern):
                        has_piperazine_nucleophile = True
                        print("Found piperazine nucleophile")
                except:
                    continue

            # Check for final product with piperazine-linked amide
            final_product_pattern = Chem.MolFromSmarts(
                "[NH][C](=[O])[CH2][N]1[CH2][CH2][N]([CH3])[CH2][CH2]1"
            )
            try:
                product_mol = Chem.MolFromSmiles(product_smiles)
                if product_mol and product_mol.HasSubstructMatch(final_product_pattern):
                    has_final_piperazine_product = True
                    print("Found final piperazine-linked product")
            except:
                pass

        # Traverse children
        for child in node.get("children", []):
            dfs_traverse(child)

    # Start traversal from the root
    dfs_traverse(route)

    # Strategy is present if all three conditions are met
    return (
        has_chloroacetamide_intermediate
        and has_piperazine_nucleophile
        and has_final_piperazine_product
    )
