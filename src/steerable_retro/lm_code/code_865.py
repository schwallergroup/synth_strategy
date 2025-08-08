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
    Detects if the synthetic route involves a late-stage Sonogashira coupling
    (C(sp²)-C(sp) bond formation between aryl halide and terminal alkyne).
    """
    sonogashira_found = False
    sonogashira_depth = -1

    def dfs_traverse(node):
        nonlocal sonogashira_found, sonogashira_depth

        if node["type"] == "reaction" and "metadata" in node and "rsmi" in node["metadata"]:
            rsmi = node["metadata"]["rsmi"]
            depth = (
                int(re.search(r"Depth: (\d+)", rsmi).group(1))
                if re.search(r"Depth: (\d+)", rsmi)
                else -1
            )

            # Only consider reactions in the second half of synthesis (late stage)
            if depth <= 3:  # Assuming depth 0-3 is late stage
                reactants = rsmi.split(">")[0].split(".")
                product = rsmi.split(">")[-1]

                try:
                    product_mol = Chem.MolFromSmiles(product)
                    reactant_mols = [Chem.MolFromSmiles(r) for r in reactants]

                    # Patterns for Sonogashira coupling
                    aryl_halide_pattern = Chem.MolFromSmarts("[c][Br,I,Cl]")
                    terminal_alkyne_pattern = Chem.MolFromSmarts("[C]#[C]")

                    # In retrosynthesis, check if product has C-C≡C and reactants have aryl halide and alkyne
                    aryl_alkyne_pattern = Chem.MolFromSmarts("[c][C]#[C]")

                    if (
                        product_mol
                        and aryl_alkyne_pattern
                        and product_mol.HasSubstructMatch(aryl_alkyne_pattern)
                    ):
                        # Check if reactants contain aryl halide and terminal alkyne
                        has_aryl_halide = False
                        has_terminal_alkyne = False

                        for r_mol in reactant_mols:
                            if r_mol:
                                if aryl_halide_pattern and r_mol.HasSubstructMatch(
                                    aryl_halide_pattern
                                ):
                                    has_aryl_halide = True
                                if terminal_alkyne_pattern and r_mol.HasSubstructMatch(
                                    terminal_alkyne_pattern
                                ):
                                    has_terminal_alkyne = True

                        if has_aryl_halide and has_terminal_alkyne:
                            sonogashira_found = True
                            sonogashira_depth = depth
                            print(f"Found Sonogashira coupling at depth {depth}: {rsmi}")
                except:
                    print(f"Error processing reaction SMILES: {rsmi}")

        # Process children
        for child in node.get("children", []):
            dfs_traverse(child)

    # Start traversal
    dfs_traverse(route)

    # Return True if Sonogashira coupling is found in late stage
    return sonogashira_found
