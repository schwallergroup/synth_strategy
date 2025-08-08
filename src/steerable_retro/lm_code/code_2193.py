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
    Detects if the synthesis involves Sonogashira coupling for C-C bond formation
    between an aryl halide and a terminal alkyne.
    """
    has_sonogashira = False

    def dfs_traverse(node):
        nonlocal has_sonogashira

        if node.get("type") == "reaction" and "metadata" in node and "rsmi" in node["metadata"]:
            rsmi = node["metadata"]["rsmi"]
            reactants = rsmi.split(">")[0].split(".")
            product = rsmi.split(">")[-1]

            # Create molecules
            product_mol = Chem.MolFromSmiles(product)
            reactant_mols = [Chem.MolFromSmiles(r) for r in reactants]

            if product_mol and all(reactant_mols):
                # Check for Sonogashira pattern: aryl halide + terminal alkyne
                aryl_halide_pattern = Chem.MolFromSmarts("c[Br,I,Cl]")
                terminal_alkyne_pattern = Chem.MolFromSmarts("C#C[H]")
                tms_alkyne_pattern = Chem.MolFromSmarts("C#C[Si]")

                has_aryl_halide = False
                has_terminal_alkyne = False
                has_tms_alkyne = False

                for r_mol in reactant_mols:
                    if r_mol.HasSubstructMatch(aryl_halide_pattern):
                        has_aryl_halide = True
                    if r_mol.HasSubstructMatch(terminal_alkyne_pattern):
                        has_terminal_alkyne = True
                    if r_mol.HasSubstructMatch(tms_alkyne_pattern):
                        has_tms_alkyne = True

                # Check if product has C-C≡C pattern (aryl-alkyne)
                aryl_alkyne_pattern = Chem.MolFromSmarts("c-C#C")
                has_aryl_alkyne_product = product_mol.HasSubstructMatch(aryl_alkyne_pattern)

                # If we have an aryl halide and a terminal/TMS alkyne as reactants,
                # and an aryl-alkyne product, it's likely a Sonogashira coupling
                if (
                    has_aryl_halide
                    and (has_terminal_alkyne or has_tms_alkyne)
                    and has_aryl_alkyne_product
                ):
                    has_sonogashira = True
                    print(f"Sonogashira coupling detected: {rsmi}")

        # Traverse children
        for child in node.get("children", []):
            dfs_traverse(child)

    # Start traversal
    dfs_traverse(route)

    return has_sonogashira
