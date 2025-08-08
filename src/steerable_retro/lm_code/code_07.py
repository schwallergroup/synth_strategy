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
    This function detects if the heterocyclic cores (pyrazole and pyridine) are preserved
    throughout the synthesis.
    """
    cores_preserved = True

    def dfs_traverse(node):
        nonlocal cores_preserved

        if node["type"] == "reaction" and "metadata" in node and "rsmi" in node["metadata"]:
            rsmi = node["metadata"]["rsmi"]
            reactants_part = rsmi.split(">")[0]
            product_part = rsmi.split(">")[-1]

            # Define patterns for pyrazole and pyridine cores
            pyrazole_pattern = Chem.MolFromSmarts("[n]1[n][c][c][c]1")
            pyridine_pattern = Chem.MolFromSmarts("[n]1[c][c][c][c][c]1")

            product_mol = Chem.MolFromSmiles(product_part)

            # Check if product has both cores
            if product_mol:
                has_pyrazole = product_mol.HasSubstructMatch(pyrazole_pattern)
                has_pyridine = product_mol.HasSubstructMatch(pyridine_pattern)

                # For each reactant, check if it has the cores
                for reactant in reactants_part.split("."):
                    reactant_mol = Chem.MolFromSmiles(reactant)
                    if (
                        reactant_mol and reactant_mol.GetNumAtoms() > 10
                    ):  # Only check complex reactants
                        reactant_has_pyrazole = reactant_mol.HasSubstructMatch(pyrazole_pattern)
                        reactant_has_pyridine = reactant_mol.HasSubstructMatch(pyridine_pattern)

                        # If cores are present in reactant but not in product, they're not preserved
                        if (reactant_has_pyrazole and not has_pyrazole) or (
                            reactant_has_pyridine and not has_pyridine
                        ):
                            print("Heterocyclic cores not preserved")
                            cores_preserved = False

        # Continue traversing
        for child in node.get("children", []):
            dfs_traverse(child)

    # Start traversal from the root
    dfs_traverse(route)
    return cores_preserved
