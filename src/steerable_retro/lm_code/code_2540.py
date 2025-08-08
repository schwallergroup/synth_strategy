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
    This function detects a specific sequence of functional group interconversions:
    carboxylic acid -> Weinreb amide -> ketone -> amide
    """
    acid_pattern = Chem.MolFromSmarts("[C](=[O])[OH]")
    weinreb_pattern = Chem.MolFromSmarts("[C](=[O])[N]([CH3])[O][CH3]")
    ketone_pattern = Chem.MolFromSmarts("[#6][C](=[O])[#6]")
    amide_pattern = Chem.MolFromSmarts("[C](=[O])[NH][#6]")

    # Track if each transformation is found
    acid_to_weinreb = False
    weinreb_to_ketone = False
    ketone_to_amide = False

    def dfs_traverse(node):
        nonlocal acid_to_weinreb, weinreb_to_ketone, ketone_to_amide

        if node["type"] == "reaction":
            if "rsmi" in node.get("metadata", {}):
                rsmi = node["metadata"]["rsmi"]
                reactants_smiles = rsmi.split(">")[0].split(".")
                product_smiles = rsmi.split(">")[-1]

                reactants_mols = [Chem.MolFromSmiles(r) for r in reactants_smiles if r]
                product_mol = Chem.MolFromSmiles(product_smiles) if product_smiles else None

                # Check for acid to Weinreb amide
                if product_mol and product_mol.HasSubstructMatch(weinreb_pattern):
                    has_acid = any(
                        r_mol and r_mol.HasSubstructMatch(acid_pattern)
                        for r_mol in reactants_mols
                        if r_mol
                    )
                    if has_acid:
                        acid_to_weinreb = True
                        print("Detected acid to Weinreb amide conversion")

                # Check for Weinreb amide to ketone
                if product_mol and product_mol.HasSubstructMatch(ketone_pattern):
                    has_weinreb = any(
                        r_mol and r_mol.HasSubstructMatch(weinreb_pattern)
                        for r_mol in reactants_mols
                        if r_mol
                    )
                    if has_weinreb:
                        weinreb_to_ketone = True
                        print("Detected Weinreb amide to ketone conversion")

                # Check for ketone to amide (via multi-component reaction)
                if product_mol and product_mol.HasSubstructMatch(amide_pattern):
                    has_ketone = any(
                        r_mol and r_mol.HasSubstructMatch(ketone_pattern)
                        for r_mol in reactants_mols
                        if r_mol
                    )
                    if has_ketone:
                        ketone_to_amide = True
                        print("Detected ketone to amide conversion")

        # Traverse children
        for child in node.get("children", []):
            dfs_traverse(child)

    # Start traversal
    dfs_traverse(route)

    # Return True if the complete sequence is detected
    return acid_to_weinreb and weinreb_to_ketone and ketone_to_amide
