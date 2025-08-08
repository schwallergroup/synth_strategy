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
    This function detects if the synthetic route employs a biaryl formation strategy.
    It looks for reactions forming C-C bonds between two aromatic rings.
    """
    found_biaryl = False

    def dfs_traverse(node):
        nonlocal found_biaryl

        if node["type"] == "reaction" and "metadata" in node and "rsmi" in node["metadata"]:
            rsmi = node["metadata"]["rsmi"]
            reactants = rsmi.split(">")[0].split(".")
            product = rsmi.split(">")[-1]

            # Check for aromatic rings in reactants
            aromatic_pattern = Chem.MolFromSmarts("[c]")
            # Check for biaryl bond in product
            biaryl_pattern = Chem.MolFromSmarts("[c]-[c]")

            aromatic_reactants = 0
            for reactant in reactants:
                try:
                    mol = Chem.MolFromSmiles(reactant)
                    if mol and mol.HasSubstructMatch(aromatic_pattern):
                        aromatic_reactants += 1
                except:
                    continue

            try:
                prod_mol = Chem.MolFromSmiles(product)
                if (
                    aromatic_reactants >= 2
                    and prod_mol
                    and prod_mol.HasSubstructMatch(biaryl_pattern)
                ):
                    print("Found biaryl formation reaction:", rsmi)
                    found_biaryl = True
            except:
                pass

        for child in node.get("children", []):
            dfs_traverse(child)

    dfs_traverse(route)
    return found_biaryl
