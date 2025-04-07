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
    This function detects a sequence of transformations at a benzylic position,
    specifically looking for alcohol → halide → nitrile sequence.
    """
    # Track molecules with specific functional groups at benzylic positions
    benzylic_alcohols = []
    benzylic_halides = []
    benzylic_nitriles = []

    # Track if we've found the complete sequence
    sequence_found = False

    def dfs_traverse(node):
        nonlocal sequence_found

        if node["type"] == "mol":
            try:
                smiles = node["smiles"]
                mol = Chem.MolFromSmiles(smiles)

                if mol:
                    # Check for benzylic alcohol
                    if mol.HasSubstructMatch(Chem.MolFromSmarts("[c][CH2][OH]")):
                        benzylic_alcohols.append(smiles)

                    # Check for benzylic halide
                    if mol.HasSubstructMatch(Chem.MolFromSmarts("[c][CH2][Br,Cl,I,F]")):
                        benzylic_halides.append(smiles)

                    # Check for benzylic nitrile
                    if mol.HasSubstructMatch(Chem.MolFromSmarts("[c][CH2][C]#[N]")):
                        benzylic_nitriles.append(smiles)
            except:
                print(f"Error processing molecule: {node['smiles']}")

        elif node["type"] == "reaction" and "metadata" in node and "rsmi" in node["metadata"]:
            rsmi = node["metadata"]["rsmi"]
            reactants = rsmi.split(">")[0].split(".")
            product = rsmi.split(">")[-1]

            try:
                # Check for alcohol to halide transformation
                alcohol_to_halide = False
                for reactant in reactants:
                    reactant_mol = Chem.MolFromSmiles(reactant)
                    if reactant_mol and reactant_mol.HasSubstructMatch(
                        Chem.MolFromSmarts("[c][CH2][OH]")
                    ):
                        product_mol = Chem.MolFromSmiles(product)
                        if product_mol and product_mol.HasSubstructMatch(
                            Chem.MolFromSmarts("[c][CH2][Br,Cl,I,F]")
                        ):
                            alcohol_to_halide = True
                            print(f"Benzylic alcohol to halide transformation detected: {rsmi}")

                # Check for halide to nitrile transformation
                halide_to_nitrile = False
                for reactant in reactants:
                    reactant_mol = Chem.MolFromSmiles(reactant)
                    if reactant_mol and reactant_mol.HasSubstructMatch(
                        Chem.MolFromSmarts("[c][CH2][Br,Cl,I,F]")
                    ):
                        product_mol = Chem.MolFromSmiles(product)
                        if product_mol and product_mol.HasSubstructMatch(
                            Chem.MolFromSmarts("[c][CH2][C]#[N]")
                        ):
                            halide_to_nitrile = True
                            print(f"Benzylic halide to nitrile transformation detected: {rsmi}")

                # If we've found both transformations, we have the sequence
                if alcohol_to_halide and halide_to_nitrile:
                    sequence_found = True
            except:
                print(f"Error analyzing reaction: {rsmi}")

        # Process children
        for child in node.get("children", []):
            dfs_traverse(child)

    # Start traversal
    dfs_traverse(route)

    # If we haven't found the complete sequence directly, check if we have all components
    if not sequence_found:
        sequence_found = (
            len(benzylic_alcohols) > 0 and len(benzylic_halides) > 0 and len(benzylic_nitriles) > 0
        )

    print(f"Benzylic functional group sequence detection:")
    print(f"  Benzylic alcohols found: {len(benzylic_alcohols)}")
    print(f"  Benzylic halides found: {len(benzylic_halides)}")
    print(f"  Benzylic nitriles found: {len(benzylic_nitriles)}")
    print(f"  Complete sequence detected: {sequence_found}")

    return sequence_found
