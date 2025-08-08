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
    Detects if the synthetic route involves coupling of two heterocyclic fragments,
    specifically looking for benzofuran and pyridine-containing structures.
    """
    found_pattern = False

    def dfs_traverse(node, depth=0):
        nonlocal found_pattern

        if node["type"] == "reaction":
            rsmi = node["metadata"]["rsmi"]
            reactants = rsmi.split(">")[0].split(".")
            product = rsmi.split(">")[-1]

            try:
                reactant_mols = [Chem.MolFromSmiles(r) for r in reactants]
                product_mol = Chem.MolFromSmiles(product)

                if product_mol is None or any(mol is None for mol in reactant_mols):
                    return

                # Define patterns for benzofuran and pyridine
                benzofuran_pattern = Chem.MolFromSmarts("c1cc2oc(=O)nc2cc1")  # Simplified pattern
                pyridine_pattern = Chem.MolFromSmarts("c1cncc([F])c1")  # Pyridine with fluorine

                # Check if product contains both heterocycles
                product_has_benzofuran = product_mol.HasSubstructMatch(benzofuran_pattern)
                product_has_pyridine = product_mol.HasSubstructMatch(pyridine_pattern)

                # Check if reactants separately contain the heterocycles
                reactants_have_benzofuran = any(
                    mol is not None and mol.HasSubstructMatch(benzofuran_pattern)
                    for mol in reactant_mols
                )
                reactants_have_pyridine = any(
                    mol is not None and mol.HasSubstructMatch(pyridine_pattern)
                    for mol in reactant_mols
                )

                # If product has both but individual reactants have one each, it's a heterocycle coupling
                if (
                    product_has_benzofuran
                    and product_has_pyridine
                    and reactants_have_benzofuran
                    and reactants_have_pyridine
                ):
                    print(f"Found heterocycle coupling at depth {depth}")
                    found_pattern = True
            except:
                print("Error processing reaction SMILES")

        # Continue traversal
        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    # Start traversal from the root
    dfs_traverse(route)
    return found_pattern
