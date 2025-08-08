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
    Detects if the synthesis involves heterocycle formation in early stages
    and convergent coupling in late stage.
    """
    has_heterocycle_formation = False
    has_late_coupling = False

    def dfs_traverse(node, depth=0):
        nonlocal has_heterocycle_formation, has_late_coupling

        if node["type"] == "reaction":
            # Check if this is a reaction node
            if "metadata" in node and "rsmi" in node["metadata"]:
                rsmi = node["metadata"]["rsmi"]
                reactants = rsmi.split(">")[0].split(".")
                product = rsmi.split(">")[-1]

                # Check for heterocycle formation in early stages (depth > 2)
                if depth > 2:
                    # Check if product has a pyrazole or similar heterocycle
                    pyrazole_pattern = Chem.MolFromSmarts("[#6]1[#6][#7][#7][#6]1")

                    try:
                        product_mol = Chem.MolFromSmiles(product)
                        if product_mol and product_mol.HasSubstructMatch(pyrazole_pattern):
                            # Check if reactants don't have the heterocycle
                            has_reactant_heterocycle = False
                            for reactant in reactants:
                                try:
                                    reactant_mol = Chem.MolFromSmiles(reactant)
                                    if reactant_mol and reactant_mol.HasSubstructMatch(
                                        pyrazole_pattern
                                    ):
                                        has_reactant_heterocycle = True
                                except:
                                    continue

                            if not has_reactant_heterocycle:
                                has_heterocycle_formation = True
                                print(f"Detected heterocycle formation at depth {depth}")
                    except:
                        pass

                # Check for convergent coupling in late stage (depth 0-1)
                if depth <= 1:
                    # Look for C-C bond formation between two complex fragments
                    # Simplified check: look for coupling reactions with complex reactants
                    complex_reactants = 0
                    for reactant in reactants:
                        # Count atoms as a simple complexity measure
                        atom_count = reactant.count(":")
                        if atom_count > 10:  # Arbitrary threshold for "complex"
                            complex_reactants += 1

                    if complex_reactants >= 2:
                        has_late_coupling = True
                        print(f"Detected convergent coupling at depth {depth}")

        # Traverse children
        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    # Start traversal
    dfs_traverse(route)

    # Return True if both conditions are met
    return has_heterocycle_formation and has_late_coupling
