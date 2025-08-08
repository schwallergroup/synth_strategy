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
    This function detects a strategy involving coupling of aromatic fragments through
    ether linkages.
    """
    # Track aromatic fragment couplings
    aromatic_couplings = 0

    def dfs_traverse(node):
        nonlocal aromatic_couplings

        if node["type"] == "reaction":
            # Extract reactants and products
            rsmi = node["metadata"]["rsmi"]
            reactants_smiles = rsmi.split(">")[0].split(".")
            product_smiles = rsmi.split(">")[-1]

            # Check if reaction involves connecting aromatic fragments
            reactant_mols = [Chem.MolFromSmiles(smi) for smi in reactants_smiles if smi]
            product_mol = Chem.MolFromSmiles(product_smiles) if product_smiles else None

            # Count aromatic rings in reactants and product
            if product_mol and len(reactant_mols) > 1:
                aromatic_rings_in_product = len(Chem.GetSSSR(product_mol))
                aromatic_rings_in_reactants = sum(
                    len(Chem.GetSSSR(mol)) for mol in reactant_mols if mol
                )

                # If product has same number of aromatic rings as reactants combined,
                # and contains ether linkages, it's likely an aromatic fragment coupling
                if aromatic_rings_in_product == aromatic_rings_in_reactants:
                    # Check for ether linkages in product
                    ether_pattern = Chem.MolFromSmarts("[c][O][c]")
                    if product_mol.HasSubstructMatch(ether_pattern):
                        print(f"Detected aromatic fragment coupling via ether in reaction: {rsmi}")
                        aromatic_couplings += 1

        # Traverse children
        for child in node.get("children", []):
            dfs_traverse(child)

    # Start traversal
    dfs_traverse(route)

    # Strategy is present if we have at least one aromatic fragment coupling
    strategy_present = aromatic_couplings > 0
    print(f"Aromatic fragment coupling strategy detected: {strategy_present}")
    print(f"Aromatic fragment couplings: {aromatic_couplings}")

    return strategy_present
