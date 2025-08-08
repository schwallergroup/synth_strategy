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
    This function detects a synthetic strategy involving nitro reduction to amine
    followed by thiourea formation via isothiocyanate coupling.
    """
    # Track if we've found each transformation
    found_nitro_reduction = False
    found_thiourea_formation = False

    def dfs_traverse(node):
        nonlocal found_nitro_reduction, found_thiourea_formation

        if node["type"] == "reaction":
            # Extract reactants and product from reaction SMILES
            rsmi = node["metadata"]["rsmi"]
            reactants_smiles = rsmi.split(">")[0].split(".")
            product_smiles = rsmi.split(">")[-1]

            # Check for nitro reduction (NO2 → NH2)
            product_mol = Chem.MolFromSmiles(product_smiles)
            if product_mol:
                amine_pattern = Chem.MolFromSmarts("[NH2]c")
                if product_mol.HasSubstructMatch(amine_pattern):
                    for reactant in reactants_smiles:
                        reactant_mol = Chem.MolFromSmiles(reactant)
                        if reactant_mol:
                            nitro_pattern = Chem.MolFromSmarts("[N+](=O)[O-]")
                            if reactant_mol.HasSubstructMatch(nitro_pattern):
                                found_nitro_reduction = True
                                print("Found nitro reduction")
                                break

            # Check for thiourea formation via isothiocyanate
            product_mol = Chem.MolFromSmiles(product_smiles)
            if product_mol:
                thiourea_pattern = Chem.MolFromSmarts("NC(=S)N")
                if product_mol.HasSubstructMatch(thiourea_pattern):
                    for reactant in reactants_smiles:
                        reactant_mol = Chem.MolFromSmiles(reactant)
                        if reactant_mol:
                            isothiocyanate_pattern = Chem.MolFromSmarts("N=C=S")
                            if reactant_mol.HasSubstructMatch(isothiocyanate_pattern):
                                found_thiourea_formation = True
                                print("Found thiourea formation via isothiocyanate")
                                break

        # Recursively traverse children
        for child in node.get("children", []):
            dfs_traverse(child)

    # Start traversal from the root
    dfs_traverse(route)

    # Return True if both transformations are found
    return found_nitro_reduction and found_thiourea_formation
