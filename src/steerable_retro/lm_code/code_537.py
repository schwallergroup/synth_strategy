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
    This function detects a linear synthesis strategy involving sequential ether bond
    formations and cleavages without convergent steps.
    """
    # Track ether bond operations and synthesis structure
    ether_operations = 0
    max_branching = 0

    # SMARTS patterns
    ether_pattern = Chem.MolFromSmarts("[c][O][c]")

    def dfs_traverse(node):
        nonlocal ether_operations, max_branching

        if node["type"] == "reaction":
            # Extract reactants and products
            rsmi = node["metadata"]["rsmi"]
            reactants_smiles = rsmi.split(">")[0].split(".")
            product_smiles = rsmi.split(">")[-1]

            # Count number of reactants to assess branching
            num_reactants = len([r for r in reactants_smiles if r])
            max_branching = max(max_branching, num_reactants)

            # Check for ether bond operations
            reactant_mols = [Chem.MolFromSmiles(smi) for smi in reactants_smiles if smi]
            product_mol = Chem.MolFromSmiles(product_smiles) if product_smiles else None

            if product_mol and product_mol.HasSubstructMatch(ether_pattern):
                reactants_with_ether = sum(
                    1 for mol in reactant_mols if mol and mol.HasSubstructMatch(ether_pattern)
                )
                if reactants_with_ether < len(reactant_mols):
                    print(f"Detected ether bond formation in reaction: {rsmi}")
                    ether_operations += 1

            if any(mol and mol.HasSubstructMatch(ether_pattern) for mol in reactant_mols) and (
                not product_mol or not product_mol.HasSubstructMatch(ether_pattern)
            ):
                print(f"Detected ether bond cleavage in reaction: {rsmi}")
                ether_operations += 1

        # Traverse children
        for child in node.get("children", []):
            dfs_traverse(child)

    # Start traversal
    dfs_traverse(route)

    # Strategy is present if we have ether operations and linear synthesis (max_branching <= 2)
    strategy_present = ether_operations > 0 and max_branching <= 2
    print(f"Linear ether synthesis strategy detected: {strategy_present}")
    print(f"Ether bond operations: {ether_operations}")
    print(f"Maximum branching factor: {max_branching}")

    return strategy_present
