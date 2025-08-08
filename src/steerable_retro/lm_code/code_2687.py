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
    This function detects the formation of a urea moiety, particularly N-hydroxyl urea.
    """
    urea_formed = False

    def dfs_traverse(node, depth=0):
        nonlocal urea_formed

        if node["type"] == "reaction":
            rsmi = node["metadata"].get("rsmi", "")
            if not rsmi:
                return

            reactants_smiles = rsmi.split(">")[0].split(".")
            product_smiles = rsmi.split(">")[-1]

            reactants_mols = [Chem.MolFromSmiles(r) for r in reactants_smiles if r]
            product_mol = Chem.MolFromSmiles(product_smiles) if product_smiles else None

            if reactants_mols and product_mol:
                # Check for N-hydroxyl urea formation
                urea_pattern = Chem.MolFromSmarts("[N]-C(=O)-[N]-[OH]")

                # Check if urea is in product but not in reactants
                product_has_urea = product_mol.HasSubstructMatch(urea_pattern)
                reactants_have_urea = any(
                    r and r.HasSubstructMatch(urea_pattern) for r in reactants_mols if r
                )

                if product_has_urea and not reactants_have_urea:
                    urea_formed = True
                    print(f"N-hydroxyl urea formation detected at depth {depth}")

        # Traverse children
        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    # Start traversal
    dfs_traverse(route)

    return urea_formed
