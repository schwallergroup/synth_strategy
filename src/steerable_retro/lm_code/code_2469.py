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
    This function detects a synthetic strategy involving late-stage nitro reduction.
    """
    has_nitro_reduction = False

    def dfs_traverse(node):
        nonlocal has_nitro_reduction

        if node["type"] == "reaction" and node.get("depth", 0) <= 1:  # Late stage (low depth)
            if "rsmi" in node.get("metadata", {}):
                rsmi = node["metadata"]["rsmi"]
                reactants = rsmi.split(">")[0].split(".")
                product = rsmi.split(">")[-1]

                # Check for nitro group in reactants
                nitro_pattern = Chem.MolFromSmarts("[#6]-[N+](=[O])[O-]")

                # Check for amine group in product
                amine_pattern = Chem.MolFromSmarts("[#6]-[NH2]")

                try:
                    reactant_mol = Chem.MolFromSmiles(reactants[0])
                    product_mol = Chem.MolFromSmiles(product)

                    if (
                        reactant_mol
                        and product_mol
                        and reactant_mol.HasSubstructMatch(nitro_pattern)
                        and product_mol.HasSubstructMatch(amine_pattern)
                    ):
                        has_nitro_reduction = True
                        print("Found late-stage nitro reduction")
                except:
                    pass

        for child in node.get("children", []):
            dfs_traverse(child)

    dfs_traverse(route)

    return has_nitro_reduction
