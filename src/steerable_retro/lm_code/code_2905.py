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
    Detects a strategy where an alkyne undergoes functionalization to form
    a more complex structure, particularly through addition reactions.
    """
    # Initialize tracking variables
    has_alkyne_functionalization = False
    alkyne_reactions = []

    def dfs_traverse(node):
        nonlocal has_alkyne_functionalization

        if node["type"] == "reaction":
            if "metadata" in node and "rsmi" in node["metadata"]:
                rsmi = node["metadata"]["rsmi"]
                reactants = rsmi.split(">")[0]
                product = rsmi.split(">")[-1]

                reactant_mol = Chem.MolFromSmiles(reactants)
                product_mol = Chem.MolFromSmiles(product)

                if reactant_mol and product_mol:
                    # Check for alkyne in reactants
                    if reactant_mol.HasSubstructMatch(Chem.MolFromSmarts("[C]#[C]")):
                        # Check for vinyl group in product (indicating addition)
                        if product_mol.HasSubstructMatch(Chem.MolFromSmarts("[C]=[C]")):
                            has_alkyne_functionalization = True
                            alkyne_reactions.append(rsmi)
                            print(f"Detected alkyne functionalization: {rsmi}")

        # Traverse children
        for child in node.get("children", []):
            dfs_traverse(child)

    # Start traversal
    dfs_traverse(route)

    # Check if the strategy is present
    strategy_present = has_alkyne_functionalization

    if strategy_present:
        print(f"Detected alkyne functionalization strategy with {len(alkyne_reactions)} reactions")
    else:
        print("Alkyne functionalization strategy not detected")

    return strategy_present
