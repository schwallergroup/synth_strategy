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
    Detects if the synthesis involves a protection/deprotection pattern.
    """
    protection_found = False
    deprotection_found = False

    def dfs_traverse(node):
        nonlocal protection_found, deprotection_found

        if node["type"] == "reaction" and "metadata" in node and "rsmi" in node["metadata"]:
            rsmi = node["metadata"]["rsmi"]
            reactants = rsmi.split(">")[0].split(".")
            product = rsmi.split(">")[-1]

            # Check for Boc protection
            product_mol = Chem.MolFromSmiles(product)
            if product_mol:
                boc_pattern = Chem.MolFromSmarts("[NX3]C(=O)OC(C)(C)C")
                if product_mol.HasSubstructMatch(boc_pattern):
                    for reactant in reactants:
                        reactant_mol = Chem.MolFromSmiles(reactant)
                        if reactant_mol:
                            amine_pattern = Chem.MolFromSmarts("[NX3H]")
                            if reactant_mol.HasSubstructMatch(amine_pattern):
                                protection_found = True
                                print("Found Boc protection")
                                break

            # Check for Boc deprotection
            for reactant in reactants:
                reactant_mol = Chem.MolFromSmiles(reactant)
                if reactant_mol:
                    boc_pattern = Chem.MolFromSmarts("[NX3]C(=O)OC(C)(C)C")
                    if reactant_mol.HasSubstructMatch(boc_pattern):
                        product_mol = Chem.MolFromSmiles(product)
                        if product_mol:
                            amine_pattern = Chem.MolFromSmarts("[NX3H]")
                            if product_mol.HasSubstructMatch(amine_pattern):
                                deprotection_found = True
                                print("Found Boc deprotection")
                                break

        # Traverse children
        for child in node.get("children", []):
            dfs_traverse(child)

    dfs_traverse(route)
    return protection_found or deprotection_found
