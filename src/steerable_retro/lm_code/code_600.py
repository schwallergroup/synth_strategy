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
    This function detects if the synthesis involves the formation of a benzothiophene ring
    from non-ring precursors (typically involving an aldehyde and a thiol).
    """
    # Initialize tracking variable
    benzothiophene_formed = False

    # Define SMARTS patterns
    benzothiophene_pattern = Chem.MolFromSmarts("c1ccc2sccc2c1")
    aldehyde_pattern = Chem.MolFromSmarts("[CH]=O")
    thiol_pattern = Chem.MolFromSmarts("[SH]")

    def dfs_traverse(node):
        nonlocal benzothiophene_formed

        if node["type"] == "reaction":
            if "rsmi" in node.get("metadata", {}):
                rsmi = node["metadata"]["rsmi"]
                reactants = rsmi.split(">")[0].split(".")
                product = rsmi.split(">")[-1]

                # Check for benzothiophene formation
                product_mol = Chem.MolFromSmiles(product)

                if product_mol and product_mol.HasSubstructMatch(benzothiophene_pattern):
                    # Check if reactants contain aldehyde and thiol
                    has_aldehyde = False
                    has_thiol = False

                    for reactant in reactants:
                        reactant_mol = Chem.MolFromSmiles(reactant)
                        if reactant_mol:
                            if reactant_mol.HasSubstructMatch(aldehyde_pattern):
                                has_aldehyde = True
                            if reactant_mol.HasSubstructMatch(thiol_pattern):
                                has_thiol = True

                    if has_aldehyde and has_thiol:
                        benzothiophene_formed = True
                        print("Found benzothiophene formation from aldehyde and thiol")

        # Traverse children
        for child in node.get("children", []):
            dfs_traverse(child)

    # Start traversal
    dfs_traverse(route)

    return benzothiophene_formed
