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
    This function detects a synthetic strategy involving multiple oxidation state changes
    on the same carbon atom (carboxylic acid → ester → alcohol → aldehyde → oxime → hydroxylamine → amide).
    """
    # Track oxidation states observed
    found_carboxylic_acid = False
    found_ester = False
    found_alcohol = False
    found_aldehyde = False
    found_oxime = False
    found_final_amide = False

    # Track transformations
    acid_to_ester = False
    ester_to_alcohol = False
    alcohol_to_aldehyde = False
    aldehyde_to_oxime = False

    def dfs_traverse(node):
        nonlocal found_carboxylic_acid, found_ester, found_alcohol, found_aldehyde, found_oxime, found_final_amide
        nonlocal acid_to_ester, ester_to_alcohol, alcohol_to_aldehyde, aldehyde_to_oxime

        if node["type"] == "mol":
            mol = Chem.MolFromSmiles(node["smiles"])
            if mol:
                # Check for functional groups
                if mol.HasSubstructMatch(Chem.MolFromSmarts("[c]-C(=O)-O")):
                    found_carboxylic_acid = True
                    print(f"Found carboxylic acid in molecule: {node['smiles']}")

                if mol.HasSubstructMatch(Chem.MolFromSmarts("[c]-C(=O)-OC")):
                    found_ester = True
                    print(f"Found ester in molecule: {node['smiles']}")

                if mol.HasSubstructMatch(Chem.MolFromSmarts("[c]-[CH2]-[OH]")):
                    found_alcohol = True
                    print(f"Found alcohol in molecule: {node['smiles']}")

                if mol.HasSubstructMatch(Chem.MolFromSmarts("[c]-[CH]=O")):
                    found_aldehyde = True
                    print(f"Found aldehyde in molecule: {node['smiles']}")

                if mol.HasSubstructMatch(Chem.MolFromSmarts("[c]-[CH]=[N]-[OH]")):
                    found_oxime = True
                    print(f"Found oxime in molecule: {node['smiles']}")

                if mol.HasSubstructMatch(Chem.MolFromSmarts("[C]-C(=O)-[N]-[OH]")):
                    found_final_amide = True
                    print(f"Found final N-hydroxy amide in molecule: {node['smiles']}")

        elif node["type"] == "reaction":
            if "metadata" in node and "rsmi" in node["metadata"]:
                rsmi = node["metadata"]["rsmi"]
                reactants = rsmi.split(">")[0].split(".")
                product = rsmi.split(">")[-1]

                reactant_mols = [Chem.MolFromSmiles(r) for r in reactants if r]
                product_mol = Chem.MolFromSmiles(product) if product else None

                if product_mol and any(r for r in reactant_mols):
                    # Check for specific transformations
                    if any(
                        r and r.HasSubstructMatch(Chem.MolFromSmarts("[c]-C(=O)-O"))
                        for r in reactant_mols
                        if r
                    ) and product_mol.HasSubstructMatch(Chem.MolFromSmarts("[c]-C(=O)-OC")):
                        acid_to_ester = True
                        print(f"Found acid to ester transformation: {rsmi}")

                    if any(
                        r and r.HasSubstructMatch(Chem.MolFromSmarts("[c]-C(=O)-OC"))
                        for r in reactant_mols
                        if r
                    ) and product_mol.HasSubstructMatch(Chem.MolFromSmarts("[c]-[CH2]-[OH]")):
                        ester_to_alcohol = True
                        print(f"Found ester to alcohol transformation: {rsmi}")

                    if any(
                        r and r.HasSubstructMatch(Chem.MolFromSmarts("[c]-[CH2]-[OH]"))
                        for r in reactant_mols
                        if r
                    ) and product_mol.HasSubstructMatch(Chem.MolFromSmarts("[c]-[CH]=O")):
                        alcohol_to_aldehyde = True
                        print(f"Found alcohol to aldehyde transformation: {rsmi}")

                    if any(
                        r and r.HasSubstructMatch(Chem.MolFromSmarts("[c]-[CH]=O"))
                        for r in reactant_mols
                        if r
                    ) and product_mol.HasSubstructMatch(Chem.MolFromSmarts("[c]-[CH]=[N]-[OH]")):
                        aldehyde_to_oxime = True
                        print(f"Found aldehyde to oxime transformation: {rsmi}")

        # Traverse children
        for child in node.get("children", []):
            dfs_traverse(child)

    # Start traversal
    dfs_traverse(route)

    # Check if we have at least 4 oxidation states and 3 transformations
    oxidation_states = sum(
        [
            found_carboxylic_acid,
            found_ester,
            found_alcohol,
            found_aldehyde,
            found_oxime,
            found_final_amide,
        ]
    )
    transformations = sum([acid_to_ester, ester_to_alcohol, alcohol_to_aldehyde, aldehyde_to_oxime])

    strategy_present = oxidation_states >= 4 and transformations >= 3
    print(f"Multi-oxidation state carbon strategy detected: {strategy_present}")
    print(f"Found {oxidation_states} oxidation states and {transformations} transformations")

    return strategy_present
