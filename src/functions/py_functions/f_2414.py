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
    This function detects the use of diverse linkers (ether and sulfonate ester) in the synthesis.
    """
    ether_linker_found = False
    sulfonate_linker_found = False

    def dfs_traverse(node):
        nonlocal ether_linker_found, sulfonate_linker_found

        if node["type"] == "reaction":
            # Extract reactants and product
            rsmi = node["metadata"]["rsmi"]
            reactants_smiles = rsmi.split(">")[0].split(".")
            product_smiles = rsmi.split(">")[-1]

            # Check for linker formation
            reactant_mols = [Chem.MolFromSmiles(smi) for smi in reactants_smiles]
            product_mol = Chem.MolFromSmiles(product_smiles)

            # Ether linker pattern (aryl-O-alkyl)
            ether_pattern = Chem.MolFromSmarts("c-O-[#6]-[#6]")

            # Sulfonate ester pattern
            sulfonate_pattern = Chem.MolFromSmarts("[#6]-O-S(=O)(=O)-[#6]")

            # Check if linkers are in product but not in reactants
            if product_mol:
                if product_mol.HasSubstructMatch(ether_pattern):
                    reactants_have_ether = any(
                        mol and mol.HasSubstructMatch(ether_pattern)
                        for mol in reactant_mols
                        if mol
                    )
                    if not reactants_have_ether:
                        print("Ether linker formation detected")
                        ether_linker_found = True

                if product_mol.HasSubstructMatch(sulfonate_pattern):
                    reactants_have_sulfonate = any(
                        mol and mol.HasSubstructMatch(sulfonate_pattern)
                        for mol in reactant_mols
                        if mol
                    )
                    if not reactants_have_sulfonate:
                        print("Sulfonate ester linker formation detected")
                        sulfonate_linker_found = True

        # Traverse children
        for child in node.get("children", []):
            dfs_traverse(child)

    # Start traversal
    dfs_traverse(route)
    return ether_linker_found and sulfonate_linker_found
