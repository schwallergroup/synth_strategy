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
    This function detects if the synthetic route employs multiple C-N bond formations
    (at least 3 different types: SNAr, amidation, sulfonamidation).
    """
    c_n_bond_types = set()

    def dfs_traverse(node):
        if node["type"] == "reaction" and node.get("metadata", {}).get("rsmi"):
            rsmi = node["metadata"]["rsmi"]
            reactants_smiles = rsmi.split(">")[0].split(".")
            product_smiles = rsmi.split(">")[-1]

            # Check for different types of C-N bond formations
            product_mol = Chem.MolFromSmiles(product_smiles)

            if product_mol:
                # Check for sulfonamide formation
                sulfonamide_pattern = Chem.MolFromSmarts("[#16](=[#8])(=[#8])[#7]")
                if product_mol.HasSubstructMatch(sulfonamide_pattern):
                    sulfonamide_in_reactants = False
                    for reactant in reactants_smiles:
                        reactant_mol = Chem.MolFromSmiles(reactant)
                        if reactant_mol and reactant_mol.HasSubstructMatch(
                            sulfonamide_pattern
                        ):
                            sulfonamide_in_reactants = True
                            break

                    if not sulfonamide_in_reactants:
                        c_n_bond_types.add("sulfonamidation")
                        print("Detected sulfonamide formation")

                # Check for amide formation
                amide_pattern = Chem.MolFromSmarts("[#6](=[#8])[#7]")
                if product_mol.HasSubstructMatch(amide_pattern):
                    # Check if this is a new amide formation
                    amide_count_product = len(
                        product_mol.GetSubstructMatches(amide_pattern)
                    )
                    amide_count_reactants = 0

                    for reactant in reactants_smiles:
                        reactant_mol = Chem.MolFromSmiles(reactant)
                        if reactant_mol:
                            amide_count_reactants += len(
                                reactant_mol.GetSubstructMatches(amide_pattern)
                            )

                    if amide_count_product > amide_count_reactants:
                        c_n_bond_types.add("amidation")
                        print("Detected amide formation")

                # Check for SNAr reaction (aromatic C-N bond formation with nitro group)
                snar_reactant_pattern = Chem.MolFromSmarts("[F][c][c]([N+](=[O])[O-])")
                for reactant in reactants_smiles:
                    reactant_mol = Chem.MolFromSmiles(reactant)
                    if reactant_mol and reactant_mol.HasSubstructMatch(
                        snar_reactant_pattern
                    ):
                        # Check if product has new C-N bond where F was
                        c_n_bond_types.add("SNAr")
                        print("Detected SNAr reaction")
                        break

        for child in node.get("children", []):
            dfs_traverse(child)

    dfs_traverse(route)
    return len(c_n_bond_types) >= 3
