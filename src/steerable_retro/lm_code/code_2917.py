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
    Detects a synthetic strategy where a heterocyclic scaffold (pyrazole)
    is maintained throughout the synthesis.
    """
    # Initialize tracking variables
    reactions_with_pyrazole = 0
    total_reactions = 0

    def dfs_traverse(node, depth=0):
        nonlocal reactions_with_pyrazole, total_reactions

        if node["type"] == "reaction":
            total_reactions += 1

            # Extract reactants and products
            rsmi = node["metadata"]["rsmi"]
            reactants_smiles = rsmi.split(">")[0].split(".")
            product_smiles = rsmi.split(">")[-1]

            # Convert to RDKit molecules
            reactants = [Chem.MolFromSmiles(smi) for smi in reactants_smiles if smi]
            product = Chem.MolFromSmiles(product_smiles) if product_smiles else None

            if product:
                # Check for pyrazole scaffold
                pyrazole_pattern = Chem.MolFromSmarts("c1nn(C)cc1")
                if product.HasSubstructMatch(pyrazole_pattern):
                    reactions_with_pyrazole += 1
                    print(f"Detected pyrazole scaffold at depth {depth}")

        # Recursively process children
        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    # Start traversal from the root
    dfs_traverse(route)

    # Check if pyrazole scaffold is maintained throughout
    scaffold_maintained = (total_reactions > 0) and (reactions_with_pyrazole == total_reactions)

    if scaffold_maintained:
        print("Maintained heterocyclic scaffold strategy detected")
    return scaffold_maintained
