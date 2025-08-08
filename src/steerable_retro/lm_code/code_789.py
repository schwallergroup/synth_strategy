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
    This function detects a strategy involving the formation of an ether linkage (C-O-C)
    between a phenol and a hydroxyethyl-containing fragment.
    """
    found_ether_formation = False

    def dfs_traverse(node, depth=0):
        nonlocal found_ether_formation

        if node["type"] == "reaction":
            rsmi = node["metadata"].get("rsmi", "")
            if not rsmi:
                return

            parts = rsmi.split(">")
            if len(parts) < 3:
                return

            reactants_smiles = parts[0].split(".")
            product_smiles = parts[2]

            # Check if product contains ether linkage
            product_mol = Chem.MolFromSmiles(product_smiles)
            if product_mol and product_mol.HasSubstructMatch(Chem.MolFromSmarts("c-[OX2]-[CX4]")):
                # Check if one reactant has phenol
                has_phenol = False
                has_alcohol = False

                for reactant in reactants_smiles:
                    reactant_mol = Chem.MolFromSmiles(reactant)
                    if not reactant_mol:
                        continue

                    if reactant_mol.HasSubstructMatch(Chem.MolFromSmarts("c[OX2H]")):
                        has_phenol = True
                    if reactant_mol.HasSubstructMatch(Chem.MolFromSmarts("[CX4][OX2H]")):
                        has_alcohol = True

                if has_phenol or has_alcohol:
                    print(f"Found ether linkage formation at depth {depth}")
                    found_ether_formation = True

        # Traverse children
        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    # Start traversal from root
    dfs_traverse(route)

    return found_ether_formation
