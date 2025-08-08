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
    This function detects if the synthetic route employs a mesylation strategy.
    It looks for reactions forming mesylates from alcohols.
    """
    found_mesylation = False

    def dfs_traverse(node):
        nonlocal found_mesylation

        if node["type"] == "reaction" and "metadata" in node and "rsmi" in node["metadata"]:
            rsmi = node["metadata"]["rsmi"]
            reactants = rsmi.split(">")[0].split(".")
            product = rsmi.split(">")[-1]

            # Check for alcohol pattern in reactants
            alcohol_pattern = Chem.MolFromSmarts("[C]-[OH]")
            # Check for mesylate pattern in product
            mesylate_pattern = Chem.MolFromSmarts("[C]-[O]-[S](=[O])(=[O])-[C]")

            has_alcohol = False
            for reactant in reactants:
                try:
                    mol = Chem.MolFromSmiles(reactant)
                    if mol and mol.HasSubstructMatch(alcohol_pattern):
                        has_alcohol = True
                        break
                except:
                    continue

            try:
                prod_mol = Chem.MolFromSmiles(product)
                if has_alcohol and prod_mol and prod_mol.HasSubstructMatch(mesylate_pattern):
                    print("Found mesylation reaction:", rsmi)
                    found_mesylation = True
            except:
                pass

        for child in node.get("children", []):
            dfs_traverse(child)

    dfs_traverse(route)
    return found_mesylation
