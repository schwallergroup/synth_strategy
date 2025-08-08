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
    This function detects N-benzylation of a heterocycle nitrogen.
    """
    benzyl_pattern = Chem.MolFromSmarts("[#6]c1ccccc1")
    n_benzyl_heterocycle_pattern = Chem.MolFromSmarts("[#7;R]([#6]c1ccccc1)[#6;R]")

    found_n_benzylation = False

    def dfs_traverse(node):
        nonlocal found_n_benzylation

        if node["type"] == "reaction":
            if "rsmi" in node.get("metadata", {}):
                rsmi = node["metadata"]["rsmi"]
                reactants = rsmi.split(">")[0].split(".")
                product = rsmi.split(">")[-1]

                reactant_mols = [Chem.MolFromSmiles(r) for r in reactants if r]
                product_mol = Chem.MolFromSmiles(product) if product else None

                if product_mol and reactant_mols:
                    # Check if product has N-benzyl heterocycle pattern
                    if product_mol.HasSubstructMatch(n_benzyl_heterocycle_pattern):
                        # Check if any reactant has benzyl group
                        has_benzyl = any(
                            reactant_mol.HasSubstructMatch(benzyl_pattern)
                            for reactant_mol in reactant_mols
                            if reactant_mol
                        )

                        if has_benzyl:
                            print("Found N-benzylation of heterocycle")
                            found_n_benzylation = True

        for child in node.get("children", []):
            dfs_traverse(child)

    dfs_traverse(route)
    return found_n_benzylation
