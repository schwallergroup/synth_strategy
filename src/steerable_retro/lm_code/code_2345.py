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
    This function detects an early-stage SNAr to introduce cyclopropylamine to aromatic core
    """
    # Define SMARTS patterns
    cyclopropylamine_pattern = Chem.MolFromSmarts("[NH][C]1[CH2][CH2]1")
    halogen_aromatic_pattern = Chem.MolFromSmarts("[c][Cl,Br,I,F]")

    # Track if we find the patterns and at what depth
    found_cyclopropylamine_depth = None
    found_halogen_aromatic_depth = None

    def dfs_traverse(node, depth=0):
        nonlocal found_cyclopropylamine_depth, found_halogen_aromatic_depth

        if node["type"] == "mol":
            mol = Chem.MolFromSmiles(node["smiles"])
            if mol:
                if mol.HasSubstructMatch(cyclopropylamine_pattern):
                    found_cyclopropylamine_depth = depth
                if mol.HasSubstructMatch(halogen_aromatic_pattern):
                    found_halogen_aromatic_depth = depth

        # Check for SNAr reaction pattern in reaction nodes
        if node["type"] == "reaction" and "metadata" in node and "rsmi" in node["metadata"]:
            rsmi = node["metadata"]["rsmi"]
            reactants = rsmi.split(">")[0].split(".")
            product = rsmi.split(">")[-1]

            # Check if reactants contain halogen and cyclopropylamine, and product has them connected
            reactant_mols = [Chem.MolFromSmiles(r) for r in reactants if r]
            product_mol = Chem.MolFromSmiles(product) if product else None

            if (
                product_mol
                and any(m and m.HasSubstructMatch(halogen_aromatic_pattern) for m in reactant_mols)
                and any(m and m.HasSubstructMatch(cyclopropylamine_pattern) for m in reactant_mols)
                and product_mol.HasSubstructMatch(cyclopropylamine_pattern)
            ):
                found_halogen_aromatic_depth = depth

        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    # Start traversal
    dfs_traverse(route)

    # Check if we found cyclopropylamine introduction via SNAr early in synthesis (high depth)
    if found_cyclopropylamine_depth is not None and found_halogen_aromatic_depth is not None:
        if found_cyclopropylamine_depth >= 4:  # Early stage (high depth in retrosynthesis)
            print(
                f"Found early-stage SNAr with cyclopropylamine at depth {found_cyclopropylamine_depth}"
            )
            return True

    print("Did not find evidence of early-stage SNAr with cyclopropylamine strategy")
    return False
