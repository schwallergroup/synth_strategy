#!/bin/python

"""LM-defined function for strategy description."""

from rdkit.Chem import AllChem, rdFMCS
import copy
from collections import deque
import rdkit.Chem as Chem
from rdkit.Chem import rdMolDescriptors
from rdkit.Chem import rdChemReactions
from rdkit.Chem import AllChem
from rdkit.Chem import rdFMCS
import rdkit.Chem.rdFMCS
from rdkit.Chem import AllChem, Descriptors, rdMolDescriptors
from rdkit import Chem
from rdkit.Chem import Descriptors
from rdkit.Chem import AllChem, rdMolDescriptors
from rdkit.Chem import AllChem, Descriptors, Lipinski
from rdkit.Chem import rdmolops
import re
from rdkit.Chem.Scaffolds import MurckoScaffold
from rdkit.Chem import AllChem, Descriptors
import traceback
import rdkit
from collections import Counter


def main(route):
    """
    This function detects a strategy involving sequential functional group
    interconversions at a benzylic position (carboxylic acid → alcohol →
    chloride → ether) combined with aromatic substitution.
    """
    # Track if we've found each transformation
    found_carboxylic_acid_reduction = False
    found_alcohol_to_chloride = False
    found_chloride_to_ether = False
    found_aryl_halide_to_nitrile = False

    def dfs_traverse(node):
        nonlocal found_carboxylic_acid_reduction, found_alcohol_to_chloride
        nonlocal found_chloride_to_ether, found_aryl_halide_to_nitrile

        if node["type"] == "reaction" and "metadata" in node and "rsmi" in node["metadata"]:
            rsmi = node["metadata"]["rsmi"]
            reactants = rsmi.split(">")[0].split(".")
            product = rsmi.split(">")[-1]

            # Check for carboxylic acid reduction to alcohol
            carboxylic_acid_pattern = Chem.MolFromSmarts("[c][C](=O)[OH]")
            benzylic_alcohol_pattern = Chem.MolFromSmarts("[c][C][OH]")

            # Check for alcohol to chloride conversion
            benzylic_chloride_pattern = Chem.MolFromSmarts("[c][C][Cl]")

            # Check for chloride to ether conversion
            benzylic_ether_pattern = Chem.MolFromSmarts("[c][C][O][C]")

            # Check for aryl halide to nitrile conversion
            aryl_bromide_pattern = Chem.MolFromSmarts("[c][Br]")
            aryl_nitrile_pattern = Chem.MolFromSmarts("[c][C]#[N]")

            # Process product molecule
            product_mol = Chem.MolFromSmiles(product)
            if product_mol:
                # Check reactants
                reactant_mols = [Chem.MolFromSmiles(r) for r in reactants if Chem.MolFromSmiles(r)]

                # Check for carboxylic acid reduction
                if any(
                    mol.HasSubstructMatch(carboxylic_acid_pattern) for mol in reactant_mols
                ) and product_mol.HasSubstructMatch(benzylic_alcohol_pattern):
                    found_carboxylic_acid_reduction = True
                    print("Found carboxylic acid reduction to benzylic alcohol")

                # Check for alcohol to chloride conversion
                if any(
                    mol.HasSubstructMatch(benzylic_alcohol_pattern) for mol in reactant_mols
                ) and product_mol.HasSubstructMatch(benzylic_chloride_pattern):
                    found_alcohol_to_chloride = True
                    print("Found benzylic alcohol to chloride conversion")

                # Check for chloride to ether conversion
                if any(
                    mol.HasSubstructMatch(benzylic_chloride_pattern) for mol in reactant_mols
                ) and product_mol.HasSubstructMatch(benzylic_ether_pattern):
                    found_chloride_to_ether = True
                    print("Found benzylic chloride to ether conversion")

                # Check for aryl halide to nitrile conversion
                if any(
                    mol.HasSubstructMatch(aryl_bromide_pattern) for mol in reactant_mols
                ) and product_mol.HasSubstructMatch(aryl_nitrile_pattern):
                    found_aryl_halide_to_nitrile = True
                    print("Found aryl halide to nitrile conversion")

        # Traverse children
        for child in node.get("children", []):
            dfs_traverse(child)

    # Start traversal
    dfs_traverse(route)

    # Strategy is present if we found at least 3 of the 4 transformations
    transformations_found = sum(
        [
            found_carboxylic_acid_reduction,
            found_alcohol_to_chloride,
            found_chloride_to_ether,
            found_aryl_halide_to_nitrile,
        ]
    )

    result = transformations_found >= 3
    print(f"Benzylic functional group interconversion strategy detected: {result}")
    return result
