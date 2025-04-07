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
    Detects synthesis routes involving carbamate formation.
    """
    carbamate_formation_detected = False

    def dfs_traverse(node):
        nonlocal carbamate_formation_detected

        if node["type"] == "reaction":
            rsmi = node.get("metadata", {}).get("rsmi", "")
            if not rsmi:
                return

            product_smiles = rsmi.split(">")[-1]

            try:
                product_mol = Chem.MolFromSmiles(product_smiles)

                if product_mol is not None:
                    # SMARTS pattern for carbamate group
                    carbamate_pattern = Chem.MolFromSmarts("[#7]-[#6](=[#8])-[#8]")

                    if product_mol.HasSubstructMatch(carbamate_pattern):
                        # Check if this is a new formation by looking at reactants
                        reactants_smiles = rsmi.split(">")[0].split(".")
                        reactant_mols = [
                            Chem.MolFromSmiles(r) for r in reactants_smiles if r
                        ]

                        # If any reactant has a carboxylic acid, this might be carbamate formation
                        carboxylic_acid_pattern = Chem.MolFromSmarts(
                            "[#6]-[#6](=[#8])-[#8;H1]"
                        )
                        has_acid = any(
                            mol is not None
                            and mol.HasSubstructMatch(carboxylic_acid_pattern)
                            for mol in reactant_mols
                        )

                        if has_acid:
                            carbamate_formation_detected = True
                            print("Detected carbamate formation from carboxylic acid")

            except Exception as e:
                print(f"Error in carbamate formation detection: {e}")

        # Process children
        for child in node.get("children", []):
            dfs_traverse(child)

    # Start traversal
    dfs_traverse(route)

    return carbamate_formation_detected
