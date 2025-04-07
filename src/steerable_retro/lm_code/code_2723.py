#!/bin/python

"""LM-defined function for strategy description."""

import copy
import re
from collections import deque

import rdkit
import rdkit.Chem as Chem
from rdkit import Chem
from rdkit.Chem import (
    AllChem,
    Descriptors,
    Lipinski,
    rdChemReactions,
    rdFMCS,
    rdMolDescriptors,
    rdmolops,
)
from rdkit.Chem.Scaffolds import MurckoScaffold


def main(route):
    """
    This function detects alcohol activation via mesylation (conversion to mesylate).
    """
    mesylation_found = False

    def dfs_traverse(node):
        nonlocal mesylation_found

        if node["type"] == "reaction":
            if "rsmi" in node.get("metadata", {}):
                rsmi = node["metadata"]["rsmi"]
                reactants = rsmi.split(">")[0].split(".")
                product = rsmi.split(">")[-1]

                # Check for mesylation
                alcohol_pattern = Chem.MolFromSmarts("[#6]-[#8;H1]")
                mesylate_pattern = Chem.MolFromSmarts("[#6]-[#8][#16](=[#8])(=[#8])[#6]")
                mesyl_chloride_pattern = Chem.MolFromSmarts("[Cl][#16](=[#8])(=[#8])[#6]")

                product_mol = Chem.MolFromSmiles(product)

                alcohol_found = False
                mesyl_chloride_found = False

                for reactant in reactants:
                    reactant_mol = Chem.MolFromSmiles(reactant)
                    if reactant_mol:
                        if reactant_mol.HasSubstructMatch(alcohol_pattern):
                            alcohol_found = True
                        if reactant_mol.HasSubstructMatch(mesyl_chloride_pattern):
                            mesyl_chloride_found = True

                if (
                    alcohol_found
                    and mesyl_chloride_found
                    and product_mol
                    and product_mol.HasSubstructMatch(mesylate_pattern)
                ):
                    mesylation_found = True
                    print("Alcohol mesylation detected")

        # Traverse children
        for child in node.get("children", []):
            dfs_traverse(child)

    # Start traversal
    dfs_traverse(route)

    return mesylation_found
