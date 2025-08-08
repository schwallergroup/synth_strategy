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
    This function detects a synthetic strategy involving benzyl ether formation and cleavage
    between aromatic fragments in a linear synthesis approach.
    """
    # Track if we find benzyl ether formations/cleavages
    benzyl_ether_formations = 0
    benzyl_ether_cleavages = 0

    # SMARTS patterns
    benzyl_ether_pattern = Chem.MolFromSmarts("[c][CH2][O][c]")
    phenol_pattern = Chem.MolFromSmarts("[c][OH]")
    benzyl_alcohol_pattern = Chem.MolFromSmarts("[c][CH2][OH]")

    def dfs_traverse(node):
        nonlocal benzyl_ether_formations, benzyl_ether_cleavages

        if node["type"] == "reaction":
            # Extract reactants and products
            rsmi = node["metadata"]["rsmi"]
            reactants_smiles = rsmi.split(">")[0].split(".")
            product_smiles = rsmi.split(">")[-1]

            # Check for benzyl ether formation (phenol + benzyl alcohol -> benzyl ether)
            reactant_mols = [Chem.MolFromSmiles(smi) for smi in reactants_smiles if smi]
            product_mol = Chem.MolFromSmiles(product_smiles) if product_smiles else None

            if product_mol and any(
                mol and mol.HasSubstructMatch(phenol_pattern) for mol in reactant_mols
            ):
                if product_mol.HasSubstructMatch(benzyl_ether_pattern):
                    print(f"Detected benzyl ether formation in reaction: {rsmi}")
                    benzyl_ether_formations += 1

            # Check for benzyl ether cleavage (benzyl ether -> phenol + benzyl derivative)
            if (
                any(mol and mol.HasSubstructMatch(phenol_pattern) for mol in reactant_mols)
                and any(
                    mol and mol.HasSubstructMatch(benzyl_alcohol_pattern) for mol in reactant_mols
                )
                and product_mol
                and product_mol.HasSubstructMatch(benzyl_ether_pattern)
            ):
                print(f"Detected benzyl ether cleavage in reaction: {rsmi}")
                benzyl_ether_cleavages += 1

        # Traverse children
        for child in node.get("children", []):
            dfs_traverse(child)

    # Start traversal
    dfs_traverse(route)

    # Strategy is present if we have at least one benzyl ether formation or cleavage
    strategy_present = benzyl_ether_formations > 0 or benzyl_ether_cleavages > 0
    print(f"Benzyl ether coupling strategy detected: {strategy_present}")
    print(f"Benzyl ether formations: {benzyl_ether_formations}")
    print(f"Benzyl ether cleavages: {benzyl_ether_cleavages}")

    return strategy_present
