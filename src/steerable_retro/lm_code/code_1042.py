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
    This function detects if the synthetic route involves an ester to amide transformation.
    """
    ester_pattern = Chem.MolFromSmarts("[#6][C](=[O])[O][#6]")
    amide_pattern = Chem.MolFromSmarts("[#6][C](=[O])[N]")
    transformation_detected = False

    def dfs_traverse(node):
        nonlocal transformation_detected

        if node["type"] == "reaction":
            if "metadata" in node and "rsmi" in node["metadata"]:
                rsmi = node["metadata"]["rsmi"]
                reactants_smiles = rsmi.split(">")[0]
                product_smiles = rsmi.split(">")[-1]

                # Check if reactants contain ester and product contains amide
                reactant_list = reactants_smiles.split(".")
                ester_in_reactants = False

                for r in reactant_list:
                    r_mol = Chem.MolFromSmiles(r)
                    if r_mol and r_mol.HasSubstructMatch(ester_pattern):
                        ester_in_reactants = True
                        break

                product_mol = Chem.MolFromSmiles(product_smiles)
                if (
                    ester_in_reactants
                    and product_mol
                    and product_mol.HasSubstructMatch(amide_pattern)
                ):
                    # Check if the amide is new (not present in all reactants)
                    all_reactants_have_amide = True
                    for r in reactant_list:
                        r_mol = Chem.MolFromSmiles(r)
                        if r_mol and not r_mol.HasSubstructMatch(amide_pattern):
                            all_reactants_have_amide = False
                            break

                    if not all_reactants_have_amide:
                        print("Ester to amide transformation detected")
                        transformation_detected = True

        # Traverse children
        for child in node.get("children", []):
            dfs_traverse(child)

    dfs_traverse(route)
    return transformation_detected
