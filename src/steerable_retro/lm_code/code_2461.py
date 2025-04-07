#!/bin/python

"""LM-defined function for strategy description."""

import copy
import re
from collections import deque

import rdkit
import rdkit.Chem as Chem
from rdkit import Chem
from rdkit.Chem import (
    AllChem,
    Descriptors,
    Lipinski,
    rdChemReactions,
    rdFMCS,
    rdMolDescriptors,
    rdmolops,
)
from rdkit.Chem.Scaffolds import MurckoScaffold


def main(route):
    """
    Detects TMS-alkyne protection/deprotection sequence in the synthetic route.
    """
    found_tms_protection = False
    found_tms_deprotection = False

    def dfs_traverse(node):
        nonlocal found_tms_protection, found_tms_deprotection

        if node["type"] == "reaction":
            # Extract reactants and product
            rsmi = node["metadata"]["rsmi"]
            reactants = rsmi.split(">")[0].split(".")
            product = rsmi.split(">")[-1]

            # Convert to RDKit molecules
            reactant_mols = [Chem.MolFromSmiles(r) for r in reactants if r]
            product_mol = Chem.MolFromSmiles(product) if product else None

            if not product_mol or not all(reactant_mols):
                return

            # Check for TMS-alkyne coupling
            tms_alkyne_pattern = Chem.MolFromSmarts("[C]#[C]-[Si]([C])([C])[C]")
            aryl_bromide_pattern = Chem.MolFromSmarts("[c]-[Br]")

            has_tms_alkyne = any(mol.HasSubstructMatch(tms_alkyne_pattern) for mol in reactant_mols)
            has_aryl_bromide = any(
                mol.HasSubstructMatch(aryl_bromide_pattern) for mol in reactant_mols
            )
            forms_tms_product = product_mol.HasSubstructMatch(tms_alkyne_pattern)

            if has_tms_alkyne and has_aryl_bromide and forms_tms_product:
                print("Found TMS-alkyne coupling")
                found_tms_protection = True

            # Check for TMS deprotection
            terminal_alkyne_pattern = Chem.MolFromSmarts("[C]#[C]-[#6]")

            has_tms_alkyne = any(mol.HasSubstructMatch(tms_alkyne_pattern) for mol in reactant_mols)
            forms_terminal_alkyne = product_mol.HasSubstructMatch(terminal_alkyne_pattern)

            if (
                has_tms_alkyne
                and forms_terminal_alkyne
                and not product_mol.HasSubstructMatch(tms_alkyne_pattern)
            ):
                print("Found TMS-alkyne deprotection")
                found_tms_deprotection = True

        # Continue traversing the tree
        for child in node.get("children", []):
            dfs_traverse(child)

    # Start traversal from the root
    dfs_traverse(route)

    # Check if we found both protection and deprotection
    return found_tms_protection and found_tms_deprotection
