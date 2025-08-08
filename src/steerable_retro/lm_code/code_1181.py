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
    Detects if the synthesis route includes methyl ether protection of phenol
    """
    found_protection = False

    def dfs_traverse(node):
        nonlocal found_protection

        if node["type"] == "reaction":
            if "rsmi" in node.get("metadata", {}):
                rsmi = node["metadata"]["rsmi"]
                reactants = rsmi.split(">")[0].split(".")
                product = rsmi.split(">")[-1]

                # Create RDKit mol objects
                product_mol = Chem.MolFromSmiles(product)
                reactant_mols = [Chem.MolFromSmiles(r) for r in reactants if r]

                if product_mol and all(reactant_mols):
                    # SMARTS for phenol
                    phenol_pattern = Chem.MolFromSmarts("[c][OH]")

                    # SMARTS for methyl ether
                    methyl_ether_pattern = Chem.MolFromSmarts("[c][O][CH3]")

                    # Check if reactants contain phenol
                    has_phenol = any(r.HasSubstructMatch(phenol_pattern) for r in reactant_mols)

                    # Check if product has methyl ether
                    has_methyl_ether = product_mol.HasSubstructMatch(methyl_ether_pattern)

                    if has_phenol and has_methyl_ether:
                        print("Found methyl ether protection")
                        found_protection = True

        # Continue traversal
        for child in node.get("children", []):
            dfs_traverse(child)

    # Start traversal from the root
    dfs_traverse(route)
    return found_protection
