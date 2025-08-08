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
    Detects a strategy involving SNAr reaction (nucleophilic aromatic substitution)
    """
    found_snar = False

    def dfs_traverse(node, depth=0):
        nonlocal found_snar

        if node["type"] == "reaction":
            # Extract reactants and product
            rsmi = node["metadata"]["rsmi"]
            reactants = rsmi.split(">")[0].split(".")
            product = rsmi.split(">")[-1]

            # Check for SNAr pattern: halogen leaving group on aromatic + amine nucleophile
            halogen_aromatic_patt = Chem.MolFromSmarts("[F,Cl,Br,I][c]")
            amine_patt = Chem.MolFromSmarts("[#7H,#7H2]")
            c_n_aromatic_patt = Chem.MolFromSmarts("[c][#7H]")

            # Check if reactants have the required patterns
            has_halogen_aromatic = False
            has_amine = False

            for reactant in reactants:
                react_mol = Chem.MolFromSmiles(reactant)
                if react_mol:
                    if react_mol.HasSubstructMatch(halogen_aromatic_patt):
                        has_halogen_aromatic = True
                    if react_mol.HasSubstructMatch(amine_patt):
                        has_amine = True

            # Check if product has C-N bond in aromatic system
            prod_mol = Chem.MolFromSmiles(product)
            has_c_n_aromatic = prod_mol and prod_mol.HasSubstructMatch(c_n_aromatic_patt)

            if has_halogen_aromatic and has_amine and has_c_n_aromatic:
                found_snar = True
                print(f"Found SNAr reaction at depth {depth}")

        # Traverse children
        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    dfs_traverse(route)

    return found_snar
