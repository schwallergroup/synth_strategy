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
    This function detects a synthetic strategy involving reductive amination.
    """
    has_reductive_amination = False

    def dfs_traverse(node):
        nonlocal has_reductive_amination

        if node["type"] == "reaction" and "metadata" in node and "rsmi" in node["metadata"]:
            rsmi = node["metadata"]["rsmi"]
            reactants = rsmi.split(">")[0].split(".")
            product = rsmi.split(">")[-1]

            # Check for aldehyde in reactants
            aldehyde_pattern = Chem.MolFromSmarts("[#6;H1]=[#8]")

            # Check for amine in reactants
            amine_pattern = Chem.MolFromSmarts("[#7;H1,H2]")

            # Check for C-N bond in product that wasn't in reactants
            cn_bond_pattern = Chem.MolFromSmarts("[#6][#7]")

            has_aldehyde = False
            has_amine = False
            for reactant in reactants:
                mol = Chem.MolFromSmiles(reactant)
                if mol:
                    if mol.HasSubstructMatch(aldehyde_pattern):
                        has_aldehyde = True
                    if mol.HasSubstructMatch(amine_pattern):
                        has_amine = True

            if has_aldehyde and has_amine:
                product_mol = Chem.MolFromSmiles(product)
                if product_mol and product_mol.HasSubstructMatch(cn_bond_pattern):
                    # This is a simplification - a more robust implementation would check
                    # that the C-N bond is specifically between the aldehyde carbon and the amine
                    has_reductive_amination = True
                    print(f"Detected reductive amination: {rsmi}")

        # Traverse children
        for child in node.get("children", []):
            dfs_traverse(child)

    # Start traversal
    dfs_traverse(route)

    print(f"Reductive amination strategy detected: {has_reductive_amination}")
    return has_reductive_amination
