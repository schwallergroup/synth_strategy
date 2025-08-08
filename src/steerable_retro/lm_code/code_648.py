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
    Detects a synthetic strategy where carboxylic acid protection is followed by
    fragment coupling via sulfonamide formation.
    """
    # Track if we found the key steps
    found_protection = False
    found_coupling = False

    def dfs_traverse(node):
        nonlocal found_protection, found_coupling

        if node["type"] == "reaction":
            # Extract reactants and product
            rsmi = node.get("metadata", {}).get("rsmi", "")
            if not rsmi:
                return

            reactants = rsmi.split(">")[0].split(".")
            product = rsmi.split(">")[-1]

            # Check for carboxylic acid protection (formation of ester)
            carboxylic_acid_pattern = Chem.MolFromSmarts("[C](=O)[OH]")
            ester_pattern = Chem.MolFromSmarts("[C](=O)O[C]")
            tert_butyl_ester_pattern = Chem.MolFromSmarts("[C](=O)OC(C)(C)C")

            # Check for sulfonamide formation
            sulfonyl_chloride_pattern = Chem.MolFromSmarts("[Cl][S](=O)(=O)[c]")
            sulfonamide_pattern = Chem.MolFromSmarts("[N][S](=O)(=O)[c]")

            # Convert SMILES to molecules
            product_mol = Chem.MolFromSmiles(product)
            reactant_mols = [Chem.MolFromSmiles(r) for r in reactants if r]

            if product_mol and reactant_mols:
                # Check for protection reaction
                if any(
                    mol.HasSubstructMatch(carboxylic_acid_pattern) for mol in reactant_mols
                ) and product_mol.HasSubstructMatch(ester_pattern):
                    print("Found carboxylic acid protection step")
                    found_protection = True

                # Check for sulfonamide coupling
                if any(
                    mol.HasSubstructMatch(sulfonyl_chloride_pattern) for mol in reactant_mols
                ) and product_mol.HasSubstructMatch(sulfonamide_pattern):
                    print("Found sulfonamide coupling step")
                    found_coupling = True

        # Continue traversing
        for child in node.get("children", []):
            dfs_traverse(child)

    # Start traversal
    dfs_traverse(route)

    # Return True if both protection and coupling were found
    return found_protection and found_coupling
