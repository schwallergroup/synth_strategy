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
    Detects if the synthetic route contains conversion of an alcohol to a mesylate.
    """
    alcohol_to_mesylate = False

    def dfs_traverse(node):
        nonlocal alcohol_to_mesylate

        if node["type"] == "reaction" and "metadata" in node and "rsmi" in node["metadata"]:
            rsmi = node["metadata"]["rsmi"]
            reactants = rsmi.split(">")[0].split(".")
            product = rsmi.split(">")[-1]

            # Check for alcohol to mesylate conversion (R-OH + MsCl → R-OMs)
            alcohol_pattern = Chem.MolFromSmarts("[OX2H]")
            mesyl_cl_pattern = Chem.MolFromSmarts("[Cl][SX4](=[OX1])(=[OX1])[CX4]")
            mesylate_pattern = Chem.MolFromSmarts("[OX2][SX4](=[OX1])(=[OX1])[CX4]")

            has_alcohol = any(
                Chem.MolFromSmiles(r) and Chem.MolFromSmiles(r).HasSubstructMatch(alcohol_pattern)
                for r in reactants
                if Chem.MolFromSmiles(r)
            )
            has_mesyl_cl = any(
                Chem.MolFromSmiles(r) and Chem.MolFromSmiles(r).HasSubstructMatch(mesyl_cl_pattern)
                for r in reactants
                if Chem.MolFromSmiles(r)
            )

            product_mol = Chem.MolFromSmiles(product)
            if (
                has_alcohol
                and has_mesyl_cl
                and product_mol
                and product_mol.HasSubstructMatch(mesylate_pattern)
            ):
                print("Found alcohol to mesylate conversion")
                alcohol_to_mesylate = True

        for child in node.get("children", []):
            dfs_traverse(child)

    dfs_traverse(route)
    return alcohol_to_mesylate
