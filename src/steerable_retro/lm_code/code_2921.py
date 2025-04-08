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
    This function detects if the synthesis involves an alcohol → azide → amine
    functional group interconversion sequence.
    """
    # Track if we've seen each transformation
    alcohol_to_azide = False
    azide_to_amine = False

    def dfs_traverse(node):
        nonlocal alcohol_to_azide, azide_to_amine

        if node["type"] == "reaction":
            if "rsmi" in node.get("metadata", {}):
                rsmi = node["metadata"]["rsmi"]
                reactants_smiles = rsmi.split(">")[0]
                product_smiles = rsmi.split(">")[-1]

                try:
                    # Check for alcohol to azide conversion
                    alcohol_pattern = Chem.MolFromSmarts("[#6]-[#8;H1]")
                    azide_pattern = Chem.MolFromSmarts("[#6]-[N]=[N+]=[N-]")

                    reactants_mol = Chem.MolFromSmiles(reactants_smiles)
                    product_mol = Chem.MolFromSmiles(product_smiles)

                    if reactants_mol and product_mol:
                        if (
                            reactants_mol.HasSubstructMatch(alcohol_pattern)
                            and product_mol.HasSubstructMatch(azide_pattern)
                            and not reactants_mol.HasSubstructMatch(azide_pattern)
                        ):
                            print("Detected alcohol to azide conversion")
                            alcohol_to_azide = True

                        # Check for azide to amine conversion
                        amine_pattern = Chem.MolFromSmarts("[#6]-[#7;H2]")
                        if (
                            reactants_mol.HasSubstructMatch(azide_pattern)
                            and product_mol.HasSubstructMatch(amine_pattern)
                            and not reactants_mol.HasSubstructMatch(amine_pattern)
                        ):
                            print("Detected azide to amine conversion")
                            azide_to_amine = True
                except Exception as e:
                    print(f"Error in functional group conversion analysis: {e}")

        for child in node.get("children", []):
            dfs_traverse(child)

    dfs_traverse(route)
    # Return True only if both transformations are detected
    return alcohol_to_azide and azide_to_amine
