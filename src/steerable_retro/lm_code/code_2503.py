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
    This function detects a synthetic strategy where:
    1. An aldehyde is protected early in the synthesis
    2. Complex fragments are built separately
    3. Reductive amination is used as a late-stage coupling
    """
    # Track if we found the key features
    found_aldehyde_protection = False
    found_late_reductive_amination = False

    # SMARTS patterns
    acetal_pattern = Chem.MolFromSmarts("[#6]1[#8][#6][#6][#8]1")  # 1,3-dioxolane pattern

    def is_reductive_amination(rsmi):
        """Check if a reaction is a reductive amination"""
        if not rsmi:
            return False

        reactants, products = rsmi.split(">")[0], rsmi.split(">")[-1]

        # Check for aldehyde and amine in reactants
        aldehyde_pattern = Chem.MolFromSmarts("[#6;H1]=[#8]")
        amine_pattern = Chem.MolFromSmarts("[#7;H1,H2]")

        # Check for secondary amine in product
        sec_amine_pattern = Chem.MolFromSmarts("[#6]-[#7;H0,H1]-[#6]")

        try:
            reactants_mol = Chem.MolFromSmiles(reactants)
            product_mol = Chem.MolFromSmiles(products)

            if (
                reactants_mol
                and product_mol
                and reactants_mol.HasSubstructMatch(aldehyde_pattern)
                and reactants_mol.HasSubstructMatch(amine_pattern)
                and product_mol.HasSubstructMatch(sec_amine_pattern)
            ):
                return True
        except:
            pass

        return False

    def is_aldehyde_protection(rsmi):
        """Check if a reaction is protecting an aldehyde as an acetal"""
        if not rsmi:
            return False

        reactants, products = rsmi.split(">")[0], rsmi.split(">")[-1]

        # Check for aldehyde in reactants
        aldehyde_pattern = Chem.MolFromSmarts("[#6;H1]=[#8]")

        try:
            reactants_mol = Chem.MolFromSmiles(reactants)
            product_mol = Chem.MolFromSmiles(products)

            if (
                reactants_mol
                and product_mol
                and reactants_mol.HasSubstructMatch(aldehyde_pattern)
                and product_mol.HasSubstructMatch(acetal_pattern)
            ):
                return True
        except:
            pass

        return False

    def dfs_traverse(node, depth=0):
        nonlocal found_aldehyde_protection, found_late_reductive_amination

        if node["type"] == "reaction":
            rsmi = node.get("metadata", {}).get("rsmi", "")

            # Check for reductive amination at low depth (late stage)
            if depth <= 1 and is_reductive_amination(rsmi):
                print(f"Found late-stage reductive amination at depth {depth}")
                found_late_reductive_amination = True

            # Check for aldehyde protection at high depth (early stage)
            if depth >= 1 and is_aldehyde_protection(rsmi):
                print(f"Found aldehyde protection at depth {depth}")
                found_aldehyde_protection = True

        # Traverse children
        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    # Start traversal
    dfs_traverse(route)

    # Return True if both key features were found
    return found_aldehyde_protection and found_late_reductive_amination
