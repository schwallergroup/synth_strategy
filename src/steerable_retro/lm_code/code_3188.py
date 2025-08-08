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
    This function detects a synthetic strategy that proceeds through sulfonyl chloride and
    isocyanate intermediates to form a sulfonylurea derivative.
    """
    # Initialize tracking variables
    has_sulfonyl_chloride = False
    has_sulfonamide = False
    has_isocyanate = False
    has_urea_or_sulfonylurea = False

    # SMARTS patterns for functional group detection
    sulfonyl_chloride_pattern = Chem.MolFromSmarts("[Cl][S](=[O])(=[O])[#6]")
    sulfonamide_pattern = Chem.MolFromSmarts("[NX3][S](=[O])(=[O])[#6]")
    isocyanate_pattern = Chem.MolFromSmarts("[#6][N]=[C]=[O]")
    urea_pattern = Chem.MolFromSmarts("[#7][C](=[O])[#7]")

    def dfs_traverse(node):
        nonlocal has_sulfonyl_chloride, has_sulfonamide, has_isocyanate, has_urea_or_sulfonylurea

        if node["type"] == "mol":
            # Check for functional groups in molecule nodes
            try:
                mol = Chem.MolFromSmiles(node["smiles"])
                if mol:
                    if mol.HasSubstructMatch(sulfonyl_chloride_pattern):
                        has_sulfonyl_chloride = True
                        print(f"Found sulfonyl chloride: {node['smiles']}")

                    if mol.HasSubstructMatch(sulfonamide_pattern):
                        has_sulfonamide = True
                        print(f"Found sulfonamide: {node['smiles']}")

                    if mol.HasSubstructMatch(isocyanate_pattern):
                        has_isocyanate = True
                        print(f"Found isocyanate: {node['smiles']}")

                    if mol.HasSubstructMatch(urea_pattern):
                        has_urea_or_sulfonylurea = True
                        print(f"Found urea/sulfonylurea: {node['smiles']}")
            except:
                print(f"Error processing molecule: {node['smiles']}")

        elif node["type"] == "reaction" and "metadata" in node and "rsmi" in node["metadata"]:
            # Analyze reaction
            rsmi = node["metadata"]["rsmi"]
            reactants = rsmi.split(">")[0].split(".")
            product = rsmi.split(">")[-1]

            # Check for functional groups in reactants and products
            for smiles in reactants + [product]:
                try:
                    mol = Chem.MolFromSmiles(smiles)
                    if mol:
                        if mol.HasSubstructMatch(sulfonyl_chloride_pattern):
                            has_sulfonyl_chloride = True
                            print(f"Found sulfonyl chloride in reaction: {smiles}")

                        if mol.HasSubstructMatch(sulfonamide_pattern):
                            has_sulfonamide = True
                            print(f"Found sulfonamide in reaction: {smiles}")

                        if mol.HasSubstructMatch(isocyanate_pattern):
                            has_isocyanate = True
                            print(f"Found isocyanate in reaction: {smiles}")

                        if mol.HasSubstructMatch(urea_pattern):
                            has_urea_or_sulfonylurea = True
                            print(f"Found urea/sulfonylurea in reaction: {smiles}")
                except:
                    print(f"Error processing molecule in reaction: {smiles}")

        # Traverse children
        for child in node.get("children", []):
            dfs_traverse(child)

    # Start traversal
    dfs_traverse(route)

    # Strategy is present if we found at least 3 of the 4 key functional groups
    strategy_present = (
        sum([has_sulfonyl_chloride, has_sulfonamide, has_isocyanate, has_urea_or_sulfonylurea]) >= 3
    )

    print(f"Sulfonamide-isocyanate-urea strategy detected: {strategy_present}")
    print(f"Found sulfonyl chloride: {has_sulfonyl_chloride}")
    print(f"Found sulfonamide: {has_sulfonamide}")
    print(f"Found isocyanate: {has_isocyanate}")
    print(f"Found urea/sulfonylurea: {has_urea_or_sulfonylurea}")

    return strategy_present
