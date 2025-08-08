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
    Detects if the synthesis route uses SNAr coupling for fragment connection.
    Looks for reactions where a chloro-substituted aromatic is coupled with a phenol.
    """
    chloro_aromatic_pattern = Chem.MolFromSmarts("[Cl]-[c]")
    phenol_pattern = Chem.MolFromSmarts("[#8H]-[c]")
    aryl_ether_pattern = Chem.MolFromSmarts("[c]-[#8]-[c]")
    has_snar = False

    def dfs_traverse(node):
        nonlocal has_snar

        if node["type"] == "reaction":
            # Extract reactants and product
            rsmi = node["metadata"]["rsmi"]
            product = rsmi.split(">")[-1]
            reactants = rsmi.split(">")[0].split(".")

            # Check for SNAr pattern: chloro-aromatic + phenol → aryl ether
            product_mol = Chem.MolFromSmiles(product)
            if product_mol and product_mol.HasSubstructMatch(aryl_ether_pattern):
                has_chloro = False
                has_phenol = False

                for reactant in reactants:
                    reactant_mol = Chem.MolFromSmiles(reactant)
                    if not reactant_mol:
                        continue

                    if reactant_mol.HasSubstructMatch(chloro_aromatic_pattern):
                        has_chloro = True
                    if reactant_mol.HasSubstructMatch(phenol_pattern):
                        has_phenol = True

                if has_chloro and has_phenol:
                    has_snar = True
                    print("Found SNAr coupling strategy")

        # Traverse children
        for child in node.get("children", []):
            dfs_traverse(child)

    dfs_traverse(route)
    return has_snar
