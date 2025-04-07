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
    Detects if the synthesis route includes a tosylate formation followed by
    nucleophilic substitution (alcohol to tosylate to substituted product).
    """
    found_tosylate_formation = False
    found_tosylate_substitution = False

    def dfs_traverse(node):
        nonlocal found_tosylate_formation, found_tosylate_substitution

        if (
            node["type"] == "reaction"
            and "metadata" in node
            and "rsmi" in node["metadata"]
        ):
            rsmi = node["metadata"]["rsmi"]
            reactants = rsmi.split(">")[0].split(".")
            product = rsmi.split(">")[-1]

            # Check for tosylate formation (alcohol + TsCl -> tosylate)
            alcohol_pattern = Chem.MolFromSmarts("[#6][OH]")
            tosyl_chloride_pattern = Chem.MolFromSmarts("c1ccc(S(=O)(=O)Cl)cc1")
            tosylate_pattern = Chem.MolFromSmarts("[#6]OS(=O)(=O)c1ccc(C)cc1")

            # Check for tosylate formation
            reactants_have_alcohol = False
            reactants_have_tosyl_chloride = False

            for reactant in reactants:
                try:
                    mol = Chem.MolFromSmiles(reactant)
                    if mol:
                        if mol.HasSubstructMatch(alcohol_pattern):
                            reactants_have_alcohol = True
                        if mol.HasSubstructMatch(tosyl_chloride_pattern):
                            reactants_have_tosyl_chloride = True
                except:
                    continue

            product_has_tosylate = False
            try:
                mol = Chem.MolFromSmiles(product)
                if mol and mol.HasSubstructMatch(tosylate_pattern):
                    product_has_tosylate = True
            except:
                pass

            if (
                reactants_have_alcohol
                and (reactants_have_tosyl_chloride or "Ts" in rsmi)
                and product_has_tosylate
            ):
                print("Found tosylate formation")
                found_tosylate_formation = True

            # Check for tosylate substitution
            tosylate_in_reactants = False
            for reactant in reactants:
                try:
                    mol = Chem.MolFromSmiles(reactant)
                    if mol and mol.HasSubstructMatch(tosylate_pattern):
                        tosylate_in_reactants = True
                        break
                except:
                    continue

            if tosylate_in_reactants and "N" in product and not product_has_tosylate:
                print("Found tosylate substitution")
                found_tosylate_substitution = True

        for child in node.get("children", []):
            dfs_traverse(child)

    dfs_traverse(route)
    return found_tosylate_formation and found_tosylate_substitution
