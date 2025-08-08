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
    This function detects a synthetic strategy involving biaryl coupling between
    heterocycles (specifically thiophene and pyridine rings).
    """
    has_biaryl_coupling = False

    def dfs_traverse(node):
        nonlocal has_biaryl_coupling

        if node["type"] == "reaction":
            if "metadata" in node and "rsmi" in node["metadata"]:
                rsmi = node["metadata"]["rsmi"]
                reactants = rsmi.split(">")[0].split(".")
                product = rsmi.split(">")[-1]

                # Check if this is a biaryl coupling reaction
                # Look for two separate heterocycles in reactants and a connected structure in product
                thiophene_pattern = Chem.MolFromSmarts("c1cscc1")
                pyridine_pattern = Chem.MolFromSmarts("c1ccncc1")

                # Check if reactants contain separate thiophene and pyridine
                reactant_mols = [Chem.MolFromSmiles(r) for r in reactants if r]
                has_thiophene = any(
                    mol.HasSubstructMatch(thiophene_pattern) for mol in reactant_mols if mol
                )
                has_pyridine = any(
                    mol.HasSubstructMatch(pyridine_pattern) for mol in reactant_mols if mol
                )

                # Check if product has connected thiophene-pyridine
                product_mol = Chem.MolFromSmiles(product)
                if product_mol:
                    biaryl_pattern = Chem.MolFromSmarts("c1cscc1-c1ccncc1")
                    if product_mol.HasSubstructMatch(biaryl_pattern):
                        print("Found biaryl coupling between thiophene and pyridine")
                        has_biaryl_coupling = True

                # Alternative approach: check for C-C bond formation between aromatic rings
                if has_thiophene and has_pyridine and product_mol:
                    # This is a simplified check - in a real implementation, you would need
                    # to analyze the actual bond changes between reactants and products
                    if product_mol.HasSubstructMatch(
                        thiophene_pattern
                    ) and product_mol.HasSubstructMatch(pyridine_pattern):
                        print("Found potential biaryl coupling (simplified check)")
                        has_biaryl_coupling = True

        # Traverse children
        for child in node.get("children", []):
            dfs_traverse(child)

    # Start traversal from the root
    dfs_traverse(route)
    return has_biaryl_coupling
