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
    This function detects a synthetic strategy involving tosylate formation
    followed by ether formation via tosylate displacement.
    """
    # Track if we found the required reactions
    found_tosylate_formation = False
    found_ether_formation = False

    def dfs_traverse(node):
        nonlocal found_tosylate_formation, found_ether_formation

        if node["type"] == "reaction":
            # Extract reactants and product from reaction SMILES
            rsmi = node.get("metadata", {}).get("rsmi", "")
            if not rsmi:
                return

            reactants = rsmi.split(">")[0].split(".")
            product = rsmi.split(">")[-1]

            # Check for tosylate formation (alcohol + sulfonyl chloride)
            if not found_tosylate_formation:
                alcohol_pattern = Chem.MolFromSmarts("[OH]")
                sulfonyl_chloride_pattern = Chem.MolFromSmarts("ClS(=O)(=O)c")
                tosylate_pattern = Chem.MolFromSmarts("OS(=O)(=O)c")

                has_alcohol = any(
                    Chem.MolFromSmiles(r)
                    and Chem.MolFromSmiles(r).HasSubstructMatch(alcohol_pattern)
                    for r in reactants
                    if r
                )
                has_sulfonyl_chloride = any(
                    Chem.MolFromSmiles(r)
                    and Chem.MolFromSmiles(r).HasSubstructMatch(sulfonyl_chloride_pattern)
                    for r in reactants
                    if r
                )
                product_mol = Chem.MolFromSmiles(product) if product else None
                has_tosylate_product = product_mol and product_mol.HasSubstructMatch(
                    tosylate_pattern
                )

                if has_alcohol and has_sulfonyl_chloride and has_tosylate_product:
                    found_tosylate_formation = True
                    print("Found tosylate formation reaction")

            # Check for ether formation via tosylate displacement
            if not found_ether_formation:
                phenol_pattern = Chem.MolFromSmarts("c[OH]")
                tosylate_pattern = Chem.MolFromSmarts("OS(=O)(=O)c")
                ether_pattern = Chem.MolFromSmarts("cO[#6]")

                has_phenol = any(
                    Chem.MolFromSmiles(r)
                    and Chem.MolFromSmiles(r).HasSubstructMatch(phenol_pattern)
                    for r in reactants
                    if r
                )
                has_tosylate = any(
                    Chem.MolFromSmiles(r)
                    and Chem.MolFromSmiles(r).HasSubstructMatch(tosylate_pattern)
                    for r in reactants
                    if r
                )
                product_mol = Chem.MolFromSmiles(product) if product else None
                has_ether_product = product_mol and product_mol.HasSubstructMatch(ether_pattern)

                if has_phenol and has_tosylate and has_ether_product:
                    found_ether_formation = True
                    print("Found ether formation via tosylate displacement")

        # Continue traversing the tree
        for child in node.get("children", []):
            dfs_traverse(child)

    # Start traversal from the root
    dfs_traverse(route)

    # Return True if both reactions were found
    return found_tosylate_formation and found_ether_formation
