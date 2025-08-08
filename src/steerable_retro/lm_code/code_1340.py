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
    This function detects a synthetic strategy involving sequential amide bond formations.
    """
    amide_formations = 0

    def dfs_traverse(node):
        nonlocal amide_formations

        if node["type"] == "reaction" and "metadata" in node and "rsmi" in node["metadata"]:
            rsmi = node["metadata"]["rsmi"]
            reactants = rsmi.split(">")[0].split(".")
            product = rsmi.split(">")[-1]

            # Check for amide formation pattern
            reactant_mols = [Chem.MolFromSmiles(r) for r in reactants if r]
            product_mol = Chem.MolFromSmiles(product) if product else None

            if all(mol is not None for mol in reactant_mols) and product_mol is not None:
                # Look for acid chloride in reactants
                acid_chloride_pattern = Chem.MolFromSmarts("[C](=O)[Cl]")
                amine_pattern = Chem.MolFromSmarts("[NH2]")
                amide_pattern = Chem.MolFromSmarts("[NH][C](=O)")

                has_acid_chloride = any(
                    len(mol.GetSubstructMatches(acid_chloride_pattern)) > 0 for mol in reactant_mols
                )
                has_amine = any(
                    len(mol.GetSubstructMatches(amine_pattern)) > 0 for mol in reactant_mols
                )
                has_amide_product = len(product_mol.GetSubstructMatches(amide_pattern)) > 0

                if has_acid_chloride and has_amine and has_amide_product:
                    amide_formations += 1
                    print(f"Detected amide formation: {rsmi}")

        # Traverse children
        for child in node.get("children", []):
            dfs_traverse(child)

    # Start traversal
    dfs_traverse(route)

    # Return True if multiple amide formations detected
    result = amide_formations >= 2
    print(f"Sequential amide formations detected: {result} (count: {amide_formations})")
    return result
