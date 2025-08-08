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
    This function detects a convergent synthesis strategy where two heterocyclic
    fragments are combined.
    """
    is_convergent_heterocycle = False

    def dfs_traverse(node):
        nonlocal is_convergent_heterocycle

        if node["type"] == "reaction":
            if "metadata" in node and "rsmi" in node["metadata"]:
                rsmi = node["metadata"]["rsmi"]
                reactants = rsmi.split(">")[0].split(".")
                product = rsmi.split(">")[-1]

                # Check if we have multiple reactants with heterocycles
                if len(reactants) >= 2:
                    reactant_mols = [Chem.MolFromSmiles(r) for r in reactants if r]

                    # Define patterns for different heterocycles
                    thiophene_pattern = Chem.MolFromSmarts("c1cscc1")
                    pyridine_pattern = Chem.MolFromSmarts("c1ccncc1")

                    # Count heterocycles in reactants
                    heterocycle_count = 0
                    for mol in reactant_mols:
                        if mol:
                            if mol.HasSubstructMatch(thiophene_pattern):
                                heterocycle_count += 1
                            if mol.HasSubstructMatch(pyridine_pattern):
                                heterocycle_count += 1

                    # If we have at least 2 heterocycles in reactants, it's convergent
                    if heterocycle_count >= 2:
                        # Check if product has both heterocycles connected
                        product_mol = Chem.MolFromSmiles(product)
                        if product_mol:
                            if product_mol.HasSubstructMatch(
                                thiophene_pattern
                            ) and product_mol.HasSubstructMatch(pyridine_pattern):
                                print("Found convergent heterocycle synthesis")
                                is_convergent_heterocycle = True

        # Traverse children
        for child in node.get("children", []):
            dfs_traverse(child)

    # Start traversal from the root
    dfs_traverse(route)
    return is_convergent_heterocycle
