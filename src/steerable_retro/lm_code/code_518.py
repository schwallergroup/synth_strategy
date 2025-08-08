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
    This function detects if the synthetic route involves pyrazole formation using hydrazine.
    """
    pyrazole_formed = False
    hydrazine_used = False

    def dfs_traverse(node):
        nonlocal pyrazole_formed, hydrazine_used

        if node["type"] == "reaction":
            # Extract reactants and product
            rsmi = node["metadata"]["rsmi"]
            reactants_smiles = rsmi.split(">")[0].split(".")
            product_smiles = rsmi.split(">")[-1]

            # Check for hydrazine in reactants
            for reactant in reactants_smiles:
                if reactant.strip() == "[NH2:11][NH2:12]" or "[NH2][NH2]" in reactant:
                    hydrazine_used = True
                    print(f"Found hydrazine in reaction: {rsmi}")

            # Check if product contains pyrazole
            if product_smiles:
                product_mol = Chem.MolFromSmiles(product_smiles)
                if product_mol:
                    pyrazole_pattern = Chem.MolFromSmarts("[n]1[n][c][c][c]1")
                    if product_mol.HasSubstructMatch(pyrazole_pattern):
                        # Check if reactants don't have pyrazole
                        reactants_have_pyrazole = False
                        for reactant in reactants_smiles:
                            if reactant.strip():
                                reactant_mol = Chem.MolFromSmiles(reactant)
                                if reactant_mol and reactant_mol.HasSubstructMatch(
                                    pyrazole_pattern
                                ):
                                    reactants_have_pyrazole = True

                        if not reactants_have_pyrazole:
                            pyrazole_formed = True
                            print(f"Pyrazole formation detected in reaction: {rsmi}")

        # Traverse children
        for child in node.get("children", []):
            dfs_traverse(child)

    # Start traversal from the root
    dfs_traverse(route)

    # Return True if both conditions are met
    result = pyrazole_formed and hydrazine_used
    print(f"Pyrazole formation from hydrazine strategy detected: {result}")
    return result
