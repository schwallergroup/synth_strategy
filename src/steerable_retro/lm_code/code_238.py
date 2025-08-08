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
    Detects if the synthesis involves formation of an oxazolidinone ring
    from a Boc-protected amino alcohol
    """
    oxazolidinone_formed = False

    def dfs_traverse(node):
        nonlocal oxazolidinone_formed

        if node["type"] == "reaction" and "metadata" in node and "rsmi" in node["metadata"]:
            rsmi = node["metadata"]["rsmi"]
            reactants = rsmi.split(">")[0].split(".")
            product = rsmi.split(">")[-1]

            # Check for Boc group in reactants
            boc_pattern = Chem.MolFromSmarts("[C](=[O])[O][C]([C])([C])[C]")
            # Check for amino alcohol pattern
            amino_alcohol_pattern = Chem.MolFromSmarts(
                "[NH][C](=[O])[O][C]([C])([C])[C].[OH][C][C]"
            )

            # Check for oxazolidinone pattern in product
            oxazolidinone_pattern = Chem.MolFromSmarts("[C]1[O][C](=[O])[NH][C]1")

            has_boc = False
            for reactant in reactants:
                try:
                    mol = Chem.MolFromSmiles(reactant)
                    if mol and mol.HasSubstructMatch(boc_pattern):
                        has_boc = True
                        break
                except:
                    continue

            try:
                prod_mol = Chem.MolFromSmiles(product)
                if has_boc and prod_mol and prod_mol.HasSubstructMatch(oxazolidinone_pattern):
                    print("Oxazolidinone ring formation detected")
                    oxazolidinone_formed = True
            except:
                pass

        # Continue traversing
        for child in node.get("children", []):
            dfs_traverse(child)

    # Start traversal from the root
    dfs_traverse(route)
    return oxazolidinone_formed
