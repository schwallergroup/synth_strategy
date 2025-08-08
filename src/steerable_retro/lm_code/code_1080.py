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
    This function detects if the synthesis route involves fragment coupling via aromatic amination
    (replacement of aromatic chloride with amine).
    """
    amination_detected = False

    def dfs_traverse(node):
        nonlocal amination_detected

        if node["type"] == "reaction" and "metadata" in node and "rsmi" in node["metadata"]:
            rsmi = node["metadata"]["rsmi"]
            reactants = rsmi.split(">")[0].split(".")
            product = rsmi.split(">")[-1]

            # Check if we have at least 2 reactants (indicating coupling)
            if len(reactants) >= 2:
                # Check for aromatic chloride in one reactant
                ar_cl_pattern = Chem.MolFromSmarts("[c][Cl]")
                # Check for amine in another reactant
                amine_pattern = Chem.MolFromSmarts("[N;H1,H2]")
                # Check for aromatic amine in product
                ar_amine_pattern = Chem.MolFromSmarts("[c][N]")

                has_ar_cl = False
                has_amine = False

                for reactant in reactants:
                    if reactant and Chem.MolFromSmiles(reactant):
                        reactant_mol = Chem.MolFromSmiles(reactant)
                        if reactant_mol:
                            if ar_cl_pattern and reactant_mol.HasSubstructMatch(ar_cl_pattern):
                                has_ar_cl = True
                            if amine_pattern and reactant_mol.HasSubstructMatch(amine_pattern):
                                has_amine = True

                if has_ar_cl and has_amine:
                    if product and Chem.MolFromSmiles(product) and ar_amine_pattern:
                        product_mol = Chem.MolFromSmiles(product)
                        if product_mol and product_mol.HasSubstructMatch(ar_amine_pattern):
                            amination_detected = True
                            print("Detected aromatic amination for fragment coupling")

        for child in node.get("children", []):
            dfs_traverse(child)

    dfs_traverse(route)
    return amination_detected
