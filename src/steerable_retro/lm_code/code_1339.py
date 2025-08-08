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
    This function detects a synthetic strategy involving multiple nitro group reductions.
    """
    nitro_reductions = 0

    def dfs_traverse(node):
        nonlocal nitro_reductions

        if node["type"] == "reaction" and "metadata" in node and "rsmi" in node["metadata"]:
            rsmi = node["metadata"]["rsmi"]
            reactants = rsmi.split(">")[0].split(".")
            product = rsmi.split(">")[-1]

            # Check for nitro reduction pattern
            reactant_mols = [Chem.MolFromSmiles(r) for r in reactants if r]
            product_mol = Chem.MolFromSmiles(product) if product else None

            if all(mol is not None for mol in reactant_mols) and product_mol is not None:
                nitro_pattern = Chem.MolFromSmarts("[N+](=O)[O-]")
                amine_pattern = Chem.MolFromSmarts("[NH2]")

                # Count nitro groups in reactants
                nitro_count_reactants = sum(
                    len(mol.GetSubstructMatches(nitro_pattern)) for mol in reactant_mols
                )

                # Count nitro groups in product
                nitro_count_product = len(product_mol.GetSubstructMatches(nitro_pattern))

                # Count amine groups in product
                amine_count_product = len(product_mol.GetSubstructMatches(amine_pattern))

                # If nitro groups decreased and amine groups present, likely a reduction
                if nitro_count_reactants > nitro_count_product and amine_count_product > 0:
                    nitro_reductions += 1
                    print(f"Detected nitro reduction: {rsmi}")

        # Traverse children
        for child in node.get("children", []):
            dfs_traverse(child)

    # Start traversal
    dfs_traverse(route)

    # Return True if multiple nitro reductions detected
    result = nitro_reductions >= 2
    print(f"Multiple nitro reductions detected: {result} (count: {nitro_reductions})")
    return result
