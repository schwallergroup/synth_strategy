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
    This function detects if a nitrile functional group is preserved throughout the synthesis.
    """
    nitrile_present_in_product = False
    nitrile_preserved = True

    def dfs_traverse(node, is_product=True):
        nonlocal nitrile_present_in_product, nitrile_preserved

        if node["type"] == "mol" and "smiles" in node:
            try:
                mol = Chem.MolFromSmiles(node["smiles"])
                nitrile_pattern = Chem.MolFromSmarts("[C]#[N]")
                has_nitrile = mol and mol.HasSubstructMatch(nitrile_pattern)

                if is_product:
                    nitrile_present_in_product = has_nitrile
                    print(f"Product {'has' if has_nitrile else 'does not have'} nitrile group")
            except:
                pass

        elif node["type"] == "reaction":
            if "metadata" in node and "rsmi" in node["metadata"]:
                rsmi = node["metadata"]["rsmi"]
                reactants = rsmi.split(">")[0].split(".")
                product = rsmi.split(">")[-1]

                try:
                    p_mol = Chem.MolFromSmiles(product)
                    p_has_nitrile = p_mol and p_mol.HasSubstructMatch(Chem.MolFromSmarts("[C]#[N]"))

                    r_has_nitrile = False
                    for reactant in reactants:
                        r_mol = Chem.MolFromSmiles(reactant)
                        if r_mol and r_mol.HasSubstructMatch(Chem.MolFromSmarts("[C]#[N]")):
                            r_has_nitrile = True
                            break

                    # If reactant had nitrile but product doesn't, nitrile wasn't preserved
                    if r_has_nitrile and not p_has_nitrile:
                        nitrile_preserved = False
                        print("Nitrile group was not preserved in a reaction")
                except:
                    pass

        for child in node.get("children", []):
            dfs_traverse(child, False)

    # Start traversal
    dfs_traverse(route)

    result = nitrile_present_in_product and nitrile_preserved
    print(f"Nitrile-preserving strategy detected: {result}")
    return result
