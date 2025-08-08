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
    This function detects a synthetic strategy involving N-alkylation of a heterocycle
    (specifically pyrazole).
    """
    has_heterocycle_alkylation = False

    def dfs_traverse(node):
        nonlocal has_heterocycle_alkylation

        if node["type"] == "reaction" and "metadata" in node and "rsmi" in node["metadata"]:
            rsmi = node["metadata"]["rsmi"]
            reactants = rsmi.split(">")[0].split(".")
            product = rsmi.split(">")[-1]

            try:
                product_mol = Chem.MolFromSmiles(product)

                # Check for N-alkylated pyrazole in product
                n_alkylated_pyrazole = Chem.MolFromSmarts("[n]1[n]([CH2][c])[c][c][c]1")

                # Check if any reactant has an unalkylated pyrazole
                has_unalkylated_pyrazole = False
                for reactant in reactants:
                    reactant_mol = Chem.MolFromSmiles(reactant)
                    if reactant_mol:
                        unalkylated_pyrazole = Chem.MolFromSmarts("[n]1[nH][c][c][c]1")
                        if reactant_mol.HasSubstructMatch(unalkylated_pyrazole):
                            has_unalkylated_pyrazole = True
                            break

                # Check if product has N-alkylated pyrazole
                if (
                    product_mol
                    and product_mol.HasSubstructMatch(n_alkylated_pyrazole)
                    and has_unalkylated_pyrazole
                ):
                    has_heterocycle_alkylation = True
                    print(f"Detected heterocycle alkylation in reaction: {rsmi}")
            except:
                pass

        # Process children
        for child in node.get("children", []):
            dfs_traverse(child)

    # Start traversal
    dfs_traverse(route)

    print(f"Heterocycle alkylation strategy detected: {has_heterocycle_alkylation}")
    return has_heterocycle_alkylation
