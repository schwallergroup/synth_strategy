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
    Detects late-stage diversification strategy using sequential N-alkylation
    and Suzuki coupling in the final steps of the synthesis.
    """
    # Track if we found each step in the sequence
    found_n_alkylation = False
    found_suzuki_coupling = False
    n_alkylation_depth = -1
    suzuki_depth = -1

    def dfs_traverse(node, depth=0):
        nonlocal found_n_alkylation, found_suzuki_coupling, n_alkylation_depth, suzuki_depth

        if node["type"] == "reaction":
            rsmi = node["metadata"]["rsmi"]
            reactants = rsmi.split(">")[0].split(".")
            product = rsmi.split(">")[-1]

            # Check for N-alkylation (benzyl-N bond formation)
            reactant_mols = [Chem.MolFromSmiles(r) for r in reactants]
            product_mol = Chem.MolFromSmiles(product)

            # Look for benzyl bromide pattern in reactants
            benzyl_br_pattern = Chem.MolFromSmarts("[c]-[CX4]-[Br]")
            # Look for tertiary amine with benzyl group in product
            benzyl_amine_pattern = Chem.MolFromSmarts("[NX3](-[CX4]-[c])")

            has_benzyl_br = any(
                mol is not None and mol.HasSubstructMatch(benzyl_br_pattern)
                for mol in reactant_mols
            )
            has_benzyl_amine = product_mol is not None and product_mol.HasSubstructMatch(
                benzyl_amine_pattern
            )

            if has_benzyl_br and has_benzyl_amine and not found_n_alkylation:
                found_n_alkylation = True
                n_alkylation_depth = depth
                print(f"Found N-alkylation step at depth {depth}")

            # Check for Suzuki coupling (biaryl formation)
            # Look for aryl bromide and boronic acid in reactants
            aryl_br_pattern = Chem.MolFromSmarts("[c]-[Br]")
            boronic_pattern = Chem.MolFromSmarts("[c]-[B]([O])([O])")
            # Look for biaryl in product
            biaryl_pattern = Chem.MolFromSmarts("[c]-[c]")

            has_aryl_br = any(
                mol is not None and mol.HasSubstructMatch(aryl_br_pattern) for mol in reactant_mols
            )
            has_boronic = any(
                mol is not None and mol.HasSubstructMatch(boronic_pattern) for mol in reactant_mols
            )
            has_biaryl = product_mol is not None and product_mol.HasSubstructMatch(biaryl_pattern)

            if has_aryl_br and has_boronic and has_biaryl and not found_suzuki_coupling:
                found_suzuki_coupling = True
                suzuki_depth = depth
                print(f"Found Suzuki coupling step at depth {depth}")

        # Continue traversing
        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    # Start traversal
    dfs_traverse(route)

    # Return True if both steps were found and they are in the correct order
    # and they occur in the late stage (low depth values)
    if found_n_alkylation and found_suzuki_coupling:
        # Check if N-alkylation precedes Suzuki coupling in the synthetic direction
        # (which means higher depth in retrosynthetic direction)
        if n_alkylation_depth > suzuki_depth:
            # Check if these are late-stage reactions (depth < 3)
            if suzuki_depth <= 2 and n_alkylation_depth <= 3:
                print("Confirmed late-stage diversification strategy")
                return True

    return False
