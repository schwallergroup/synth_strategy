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
    This function detects if the synthetic route involves activation of a primary
    alcohol by conversion to a halide (particularly bromide).
    """
    alcohol_halogenation_found = False

    def dfs_traverse(node):
        nonlocal alcohol_halogenation_found

        if node["type"] == "reaction" and "metadata" in node and "rsmi" in node["metadata"]:
            rsmi = node["metadata"]["rsmi"]
            reactants = rsmi.split(">")[0].split(".")
            product = rsmi.split(">")[-1]

            # Check for alcohol to halide transformation
            alcohol_pattern = Chem.MolFromSmarts("[OH][CH2][c]")
            halide_pattern = Chem.MolFromSmarts("[Br,Cl,I][CH2][c]")

            reactant_mols = [Chem.MolFromSmiles(r) for r in reactants if r]
            product_mol = Chem.MolFromSmiles(product) if product else None

            if product_mol and product_mol.HasSubstructMatch(halide_pattern):
                # Check if any reactant has alcohol
                for r_mol in reactant_mols:
                    if r_mol and r_mol.HasSubstructMatch(alcohol_pattern):
                        alcohol_halogenation_found = True
                        print(f"Found alcohol to halide conversion: {rsmi}")
                        break

        # Traverse children
        for child in node.get("children", []):
            dfs_traverse(child)

    # Start traversal from root
    dfs_traverse(route)

    print(f"Alcohol activation by halogenation strategy detected: {alcohol_halogenation_found}")
    return alcohol_halogenation_found
