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
    This function detects a synthetic strategy involving ester formation and cleavage
    as a protection/deprotection approach.
    """
    ester_formation = False
    ester_cleavage = False

    def dfs_traverse(node):
        nonlocal ester_formation, ester_cleavage

        if node["type"] == "reaction":
            if "metadata" in node and "rsmi" in node["metadata"]:
                rsmi = node["metadata"]["rsmi"]
                reactants = rsmi.split(">")[0].split(".")
                product = rsmi.split(">")[-1]

                # Check for ester formation (phenol + carboxylic acid derivative)
                phenol_pattern = Chem.MolFromSmarts("[#6]1:[#6]:[#6]:[#6]:[#6]:[#6]:1-[#8]")
                carboxylic_pattern = Chem.MolFromSmarts("[#6]-[#6](=[#8])-[#8]")
                ester_pattern = Chem.MolFromSmarts(
                    "[#6]1:[#6]:[#6]:[#6]:[#6]:[#6]:1-[#8]-[#6]-[#6](=[#8])"
                )

                # Check reactants for phenol and product for ester
                reactant_mols = [Chem.MolFromSmiles(r) for r in reactants if r]
                product_mol = Chem.MolFromSmiles(product) if product else None

                if product_mol and any(
                    mol and mol.HasSubstructMatch(phenol_pattern) for mol in reactant_mols
                ):
                    if product_mol.HasSubstructMatch(ester_pattern):
                        print("Detected ester formation from phenol")
                        ester_formation = True

                # Check for ester cleavage
                if product_mol and any(
                    mol and mol.HasSubstructMatch(ester_pattern) for mol in reactant_mols
                ):
                    if product_mol.HasSubstructMatch(
                        phenol_pattern
                    ) or product_mol.HasSubstructMatch(carboxylic_pattern):
                        print("Detected ester cleavage")
                        ester_cleavage = True

        # Traverse children
        for child in node.get("children", []):
            dfs_traverse(child)

    # Start traversal
    dfs_traverse(route)

    # Return True if both ester formation and cleavage are detected
    return ester_formation and ester_cleavage
