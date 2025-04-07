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
    This function detects the use of Weinreb amide in a synthetic route,
    specifically looking for formation and subsequent reaction to form a ketone.
    """
    weinreb_amide_formed = False
    weinreb_amide_reacted = False

    def dfs_traverse(node):
        nonlocal weinreb_amide_formed, weinreb_amide_reacted

        if node["type"] == "reaction":
            if "rsmi" in node.get("metadata", {}):
                rsmi = node["metadata"]["rsmi"]
                reactants = rsmi.split(">")[0].split(".")
                product = rsmi.split(">")[-1]

                # Check for Weinreb amide formation
                product_mol = Chem.MolFromSmiles(product)
                weinreb_pattern = Chem.MolFromSmarts("[C](=[O])[N]([CH3])[O][CH3]")
                if product_mol and product_mol.HasSubstructMatch(weinreb_pattern):
                    weinreb_amide_formed = True
                    print("Weinreb amide formation detected")

                # Check for Weinreb amide reaction to form ketone
                if weinreb_amide_formed:
                    reactant_has_weinreb = False
                    for reactant in reactants:
                        reactant_mol = Chem.MolFromSmiles(reactant)
                        if reactant_mol and reactant_mol.HasSubstructMatch(
                            weinreb_pattern
                        ):
                            reactant_has_weinreb = True
                            break

                    product_mol = Chem.MolFromSmiles(product)
                    ketone_pattern = Chem.MolFromSmarts("[C](=[O])[C]")
                    if (
                        reactant_has_weinreb
                        and product_mol
                        and product_mol.HasSubstructMatch(ketone_pattern)
                    ):
                        weinreb_amide_reacted = True
                        print("Weinreb amide reaction to ketone detected")

        for child in node.get("children", []):
            dfs_traverse(child)

    dfs_traverse(route)
    return weinreb_amide_formed and weinreb_amide_reacted
