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
    Detects if the synthetic route uses chloro displacement strategy
    for attaching heterocycles.
    """
    found_pattern = False

    def dfs_traverse(node):
        nonlocal found_pattern

        if node["type"] == "reaction" and "metadata" in node and "rsmi" in node["metadata"]:
            # Extract reactants and product
            rsmi = node["metadata"]["rsmi"]
            reactants_smiles = rsmi.split(">")[0].split(".")
            product_smiles = rsmi.split(">")[-1]

            try:
                # Check for chloro-substituted aromatic in reactants
                cl_aromatic_pattern = Chem.MolFromSmarts("[c]-[Cl]")

                # Check for heterocycle in reactants
                indole_pattern = Chem.MolFromSmarts("[c]1[c][c][c]2[nH][c][c][c]2[c]1")

                # Check if one reactant has chloro group and another has heterocycle
                has_cl_aromatic = False
                has_heterocycle = False

                for reactant in reactants_smiles:
                    reactant_mol = Chem.MolFromSmiles(reactant)
                    if not reactant_mol:
                        continue

                    if reactant_mol.HasSubstructMatch(cl_aromatic_pattern):
                        has_cl_aromatic = True

                    if reactant_mol.HasSubstructMatch(indole_pattern):
                        has_heterocycle = True

                # Check if product has C-O bond where chloro was
                if has_cl_aromatic and has_heterocycle:
                    product_mol = Chem.MolFromSmiles(product_smiles)
                    co_pattern = Chem.MolFromSmarts("[c]-[O]-[c]")

                    if product_mol and product_mol.HasSubstructMatch(co_pattern):
                        print("Found chloro displacement for heterocycle attachment")
                        found_pattern = True
            except:
                pass

        # Continue traversing
        for child in node.get("children", []):
            dfs_traverse(child)

    # Start traversal
    dfs_traverse(route)
    return found_pattern
