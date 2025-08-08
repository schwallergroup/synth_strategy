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
    Detects a synthetic route that includes ring opening of a cyclic structure.
    """
    has_ring_opening = False

    def dfs_traverse(node):
        nonlocal has_ring_opening

        if node["type"] == "reaction":
            try:
                rsmi = node["metadata"]["rsmi"]
                reactants_smiles = rsmi.split(">")[0].split(".")
                product_smiles = rsmi.split(">")[-1].split(".")

                # Convert SMILES to molecules, filtering out None values
                reactants = [Chem.MolFromSmiles(smi) for smi in reactants_smiles if smi]
                products = [Chem.MolFromSmiles(smi) for smi in product_smiles if smi]

                # Filter out None molecules (failed parsing)
                reactants = [mol for mol in reactants if mol is not None]
                products = [mol for mol in products if mol is not None]

                if reactants and products:  # Ensure we have valid molecules
                    # Count rings in reactants and products
                    reactant_rings = sum(mol.GetRingInfo().NumRings() for mol in reactants)
                    product_rings = sum(mol.GetRingInfo().NumRings() for mol in products)

                    # Check if there's a decrease in ring count
                    if product_rings < reactant_rings:
                        print(f"Found ring opening: {reactant_rings} rings → {product_rings} rings")
                        has_ring_opening = True
            except Exception as e:
                print(f"Error processing reaction: {e}")

        # Traverse children
        for child in node.get("children", []):
            dfs_traverse(child)

    # Start traversal
    dfs_traverse(route)

    print(f"Ring opening strategy detected: {has_ring_opening}")
    return has_ring_opening
