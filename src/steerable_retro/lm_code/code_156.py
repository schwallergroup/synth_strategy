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
    This function detects a synthetic strategy involving late-stage amide coupling.
    """
    late_stage_amide_found = False

    def dfs_traverse(node, depth=0):
        nonlocal late_stage_amide_found

        if node["type"] == "reaction" and depth <= 1:  # Late stage = low depth
            if "rsmi" in node.get("metadata", {}):
                rsmi = node["metadata"]["rsmi"]
                reactants = rsmi.split(">")[0]
                products = rsmi.split(">")[-1]

                # Check for amine and carboxylic acid in reactants
                reactant_mol = Chem.MolFromSmiles(reactants)
                if reactant_mol:
                    has_amine = reactant_mol.HasSubstructMatch(Chem.MolFromSmarts("[NH2]"))
                    has_carboxylic = reactant_mol.HasSubstructMatch(Chem.MolFromSmarts("C(=O)[OH]"))

                    # Check for amide in products
                    product_mol = Chem.MolFromSmiles(products)
                    if product_mol and product_mol.HasSubstructMatch(
                        Chem.MolFromSmarts("[NH]C(=O)")
                    ):
                        if has_amine and has_carboxylic:
                            print(f"Detected late-stage amide coupling at depth {depth}")
                            late_stage_amide_found = True

        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    dfs_traverse(route)
    return late_stage_amide_found
