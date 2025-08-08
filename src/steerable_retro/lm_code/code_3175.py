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
    This function detects a synthetic strategy involving Sonogashira coupling
    (aryl halide + terminal alkyne → internal alkyne).
    """
    sonogashira_detected = False

    def dfs_traverse(node):
        nonlocal sonogashira_detected

        if node["type"] == "reaction" and "metadata" in node and "rsmi" in node["metadata"]:
            rsmi = node["metadata"]["rsmi"]
            reactants = rsmi.split(">")[0].split(".")
            product = rsmi.split(">")[-1]

            # Check for terminal alkyne in reactants
            terminal_alkyne_pattern = Chem.MolFromSmarts("[C]#[CH]")
            # Check for aryl halide in reactants
            aryl_halide_pattern = Chem.MolFromSmarts("[c]-[Br,Cl,I]")
            # Check for internal alkyne in product that wasn't in reactants
            internal_alkyne_pattern = Chem.MolFromSmarts("[c]-[C]#[C]")

            reactant_mols = [Chem.MolFromSmiles(r) for r in reactants if r]
            product_mol = Chem.MolFromSmiles(product) if product else None

            if (
                product_mol
                and any(
                    reactant_mol and reactant_mol.HasSubstructMatch(terminal_alkyne_pattern)
                    for reactant_mol in reactant_mols
                )
                and any(
                    reactant_mol and reactant_mol.HasSubstructMatch(aryl_halide_pattern)
                    for reactant_mol in reactant_mols
                )
                and product_mol.HasSubstructMatch(internal_alkyne_pattern)
            ):
                print("Sonogashira coupling detected in reaction:", rsmi)
                sonogashira_detected = True

        for child in node.get("children", []):
            dfs_traverse(child)

    dfs_traverse(route)
    return sonogashira_detected
