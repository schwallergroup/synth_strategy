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
    This function detects the use of ester protection strategy in the synthesis,
    looking for both methyl and tert-butyl esters.
    """
    methyl_ester_count = 0
    tbutyl_ester_count = 0
    ester_deprotection_count = 0

    def dfs_traverse(node, depth=0):
        nonlocal methyl_ester_count, tbutyl_ester_count, ester_deprotection_count

        if node["type"] == "mol":
            mol = Chem.MolFromSmiles(node["smiles"])
            if mol:
                # Check for methyl ester
                if mol.HasSubstructMatch(Chem.MolFromSmarts("[C](=[O])[O][CH3]")):
                    methyl_ester_count += 1
                    print(f"Found methyl ester in molecule at depth {depth}")

                # Check for tert-butyl ester
                if mol.HasSubstructMatch(Chem.MolFromSmarts("[C](=[O])[O]C(C)(C)C")):
                    tbutyl_ester_count += 1
                    print(f"Found tert-butyl ester in molecule at depth {depth}")

        elif node["type"] == "reaction":
            if "rsmi" in node.get("metadata", {}):
                rsmi = node["metadata"]["rsmi"]
                reactants = rsmi.split(">")[0].split(".")
                product = rsmi.split(">")[-1]

                reactant_mols = [Chem.MolFromSmiles(r) for r in reactants]
                product_mol = Chem.MolFromSmiles(product)

                # Check for ester deprotection
                has_ester_reactant = any(
                    mol
                    and (
                        mol.HasSubstructMatch(Chem.MolFromSmarts("[C](=[O])[O][CH3]"))
                        or mol.HasSubstructMatch(
                            Chem.MolFromSmarts("[C](=[O])[O]C(C)(C)C")
                        )
                    )
                    for mol in reactant_mols
                )

                has_acid_product = product_mol and product_mol.HasSubstructMatch(
                    Chem.MolFromSmarts("[C](=[O])[OH]")
                )

                if has_ester_reactant and has_acid_product:
                    ester_deprotection_count += 1
                    print(f"Found ester deprotection at depth {depth}")

        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    dfs_traverse(route)

    # Check if both types of esters are used and if there are deprotection steps
    result = (
        methyl_ester_count > 0 or tbutyl_ester_count > 0
    ) and ester_deprotection_count >= 2
    print(f"Ester protection strategy: {result}")
    print(f"Methyl ester count: {methyl_ester_count}")
    print(f"tert-Butyl ester count: {tbutyl_ester_count}")
    print(f"Ester deprotection count: {ester_deprotection_count}")

    return result
