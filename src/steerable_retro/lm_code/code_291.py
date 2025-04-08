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
    This function detects a synthetic strategy involving a diarylamine core structure
    that undergoes sequential nitrogen-containing functional group transformations
    (oxime → hydroxylamine → amide) in the late stages of synthesis.
    """
    # Track if we've found the key transformations
    found_diarylamine = False
    found_oxime = False
    found_hydroxylamine = False
    found_n_hydroxy_amide = False

    # Track the sequence of transformations
    oxime_to_hydroxylamine = False
    hydroxylamine_to_amide = False

    def dfs_traverse(node):
        nonlocal found_diarylamine, found_oxime, found_hydroxylamine, found_n_hydroxy_amide
        nonlocal oxime_to_hydroxylamine, hydroxylamine_to_amide

        if node["type"] == "mol":
            # Check for diarylamine core
            mol = Chem.MolFromSmiles(node["smiles"])
            if mol:
                diarylamine_pattern = Chem.MolFromSmarts("[c]-[NH]-[c]")
                if mol.HasSubstructMatch(diarylamine_pattern):
                    found_diarylamine = True

                # Check for specific functional groups
                oxime_pattern = Chem.MolFromSmarts("[C]=[N]-[OH]")
                if mol.HasSubstructMatch(oxime_pattern):
                    found_oxime = True
                    print(f"Found oxime in molecule: {node['smiles']}")

                hydroxylamine_pattern = Chem.MolFromSmarts("[NH]-[OH]")
                if mol.HasSubstructMatch(hydroxylamine_pattern) and not mol.HasSubstructMatch(
                    Chem.MolFromSmarts("[C](=O)-[N]-[OH]")
                ):
                    found_hydroxylamine = True
                    print(f"Found hydroxylamine in molecule: {node['smiles']}")

                n_hydroxy_amide_pattern = Chem.MolFromSmarts("[C](=O)-[N]-[OH]")
                if mol.HasSubstructMatch(n_hydroxy_amide_pattern):
                    found_n_hydroxy_amide = True
                    print(f"Found N-hydroxy amide in molecule: {node['smiles']}")

        elif node["type"] == "reaction":
            if "metadata" in node and "rsmi" in node["metadata"]:
                rsmi = node["metadata"]["rsmi"]
                reactants = rsmi.split(">")[0].split(".")
                product = rsmi.split(">")[-1]

                # Check for oxime to hydroxylamine transformation
                reactant_mols = [Chem.MolFromSmiles(r) for r in reactants if r]
                product_mol = Chem.MolFromSmiles(product) if product else None

                if product_mol and any(r for r in reactant_mols):
                    # Check if this reaction converts oxime to hydroxylamine
                    if (
                        any(
                            r and r.HasSubstructMatch(Chem.MolFromSmarts("[C]=[N]-[OH]"))
                            for r in reactant_mols
                            if r
                        )
                        and product_mol.HasSubstructMatch(Chem.MolFromSmarts("[NH]-[OH]"))
                        and not product_mol.HasSubstructMatch(
                            Chem.MolFromSmarts("[C](=O)-[N]-[OH]")
                        )
                    ):
                        oxime_to_hydroxylamine = True
                        print(f"Found oxime to hydroxylamine transformation: {rsmi}")

                    # Check if this reaction converts hydroxylamine to N-hydroxy amide
                    if any(
                        r
                        and r.HasSubstructMatch(Chem.MolFromSmarts("[NH]-[OH]"))
                        and not r.HasSubstructMatch(Chem.MolFromSmarts("[C](=O)-[N]-[OH]"))
                        for r in reactant_mols
                        if r
                    ) and product_mol.HasSubstructMatch(Chem.MolFromSmarts("[C](=O)-[N]-[OH]")):
                        hydroxylamine_to_amide = True
                        print(f"Found hydroxylamine to N-hydroxy amide transformation: {rsmi}")

        # Traverse children
        for child in node.get("children", []):
            dfs_traverse(child)

    # Start traversal
    dfs_traverse(route)

    # Check if the strategy is present
    strategy_present = (
        found_diarylamine
        and found_oxime
        and found_hydroxylamine
        and found_n_hydroxy_amide
        and oxime_to_hydroxylamine
        and hydroxylamine_to_amide
    )

    print(
        f"Diarylamine core with sequential N transformations strategy detected: {strategy_present}"
    )
    return strategy_present
