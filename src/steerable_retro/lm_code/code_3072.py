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
    Detects if the synthesis involves introduction of an alkyne side chain
    through phenolic ether formation.
    """
    has_alkyne = False
    has_phenolic_ether_formation = False

    def dfs_traverse(node, depth=0):
        nonlocal has_alkyne, has_phenolic_ether_formation

        if node["type"] == "reaction":
            if "rsmi" in node.get("metadata", {}):
                rsmi = node["metadata"]["rsmi"]
                reactants = rsmi.split(">")[0].split(".")
                product = rsmi.split(">")[-1]

                # Check for alkyne in product
                product_mol = Chem.MolFromSmiles(product)
                if product_mol:
                    alkyne_pattern = Chem.MolFromSmarts("[CX2]#[CX2]")
                    if product_mol.HasSubstructMatch(alkyne_pattern):
                        has_alkyne = True

                # Check for phenolic ether formation
                phenol_pattern = Chem.MolFromSmarts("[OH][c]")
                ether_pattern = Chem.MolFromSmarts("[O][c]")

                # Look for phenol in reactants and ether in product
                for reactant in reactants:
                    reactant_mol = Chem.MolFromSmiles(reactant)
                    if reactant_mol and reactant_mol.HasSubstructMatch(phenol_pattern):
                        product_mol = Chem.MolFromSmiles(product)
                        if product_mol and product_mol.HasSubstructMatch(ether_pattern):
                            # Check if the other reactant has an alkyne
                            for other_reactant in reactants:
                                if other_reactant != reactant:
                                    other_mol = Chem.MolFromSmiles(other_reactant)
                                    if other_mol:
                                        alkyne_pattern = Chem.MolFromSmarts("[CX2]#[CX2]")
                                        if other_mol.HasSubstructMatch(alkyne_pattern):
                                            has_phenolic_ether_formation = True
                                            print(
                                                "Detected phenolic ether formation with alkyne side chain"
                                            )

        # Continue traversing
        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    dfs_traverse(route)
    result = has_alkyne and has_phenolic_ether_formation

    if result:
        print("Detected alkyne side chain introduction via phenolic ether formation")

    return result
