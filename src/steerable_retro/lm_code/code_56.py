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
    This function detects if the synthetic route employs an N-alkylation strategy.
    It looks for reactions forming C-N bonds between an aniline and an alkyl group.
    """
    found_n_alkylation = False

    def dfs_traverse(node):
        nonlocal found_n_alkylation

        if node["type"] == "reaction" and "metadata" in node and "rsmi" in node["metadata"]:
            rsmi = node["metadata"]["rsmi"]
            reactants = rsmi.split(">")[0].split(".")
            product = rsmi.split(">")[-1]

            # Check for aniline pattern in reactants
            aniline_pattern = Chem.MolFromSmarts("[c]-[NH2]")
            # Check for N-alkylated product
            n_alkylated_pattern = Chem.MolFromSmarts("[c]-[NH]-[C]")

            has_aniline = False
            for reactant in reactants:
                try:
                    mol = Chem.MolFromSmiles(reactant)
                    if mol and mol.HasSubstructMatch(aniline_pattern):
                        has_aniline = True
                        break
                except:
                    continue

            try:
                prod_mol = Chem.MolFromSmiles(product)
                if has_aniline and prod_mol and prod_mol.HasSubstructMatch(n_alkylated_pattern):
                    print("Found N-alkylation reaction:", rsmi)
                    found_n_alkylation = True
            except:
                pass

        for child in node.get("children", []):
            dfs_traverse(child)

    dfs_traverse(route)
    return found_n_alkylation
