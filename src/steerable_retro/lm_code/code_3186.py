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
    This function detects a synthetic strategy involving Boc deprotection
    of a nitrogen-containing heterocycle.
    """
    boc_deprotection_found = False

    def dfs_traverse(node, depth=0):
        nonlocal boc_deprotection_found

        if node["type"] == "reaction":
            rsmi = node["metadata"]["rsmi"]
            reactants_smiles = rsmi.split(">")[0].split(".")
            product_smiles = rsmi.split(">")[-1]

            reactant_mols = [Chem.MolFromSmiles(r) for r in reactants_smiles if r]
            product_mol = Chem.MolFromSmiles(product_smiles) if product_smiles else None

            if product_mol and all(r is not None for r in reactant_mols):
                # Check for Boc deprotection
                boc_pattern = Chem.MolFromSmarts("[#6]C([#6])([#6])[#8]C(=O)[#7]")
                piperidine_pattern = Chem.MolFromSmarts("[#6]1[#6][#6][#6][#7][#6]1")

                if any(
                    r.HasSubstructMatch(boc_pattern) for r in reactant_mols
                ) and not product_mol.HasSubstructMatch(boc_pattern):
                    # Check if the molecule contains a piperidine ring
                    if product_mol.HasSubstructMatch(piperidine_pattern):
                        print(f"Detected Boc deprotection of piperidine at depth {depth}")
                        boc_deprotection_found = True

        # Traverse children
        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    # Start traversal from the root
    dfs_traverse(route)

    return boc_deprotection_found
