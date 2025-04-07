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
    This function detects if the synthesis involves O-alkylation with an epoxide.
    """
    found_epoxide_alkylation = False

    def dfs_traverse(node):
        nonlocal found_epoxide_alkylation

        if node["type"] == "reaction":
            if "rsmi" in node.get("metadata", {}):
                rsmi = node["metadata"]["rsmi"]
                reactants_smiles = rsmi.split(">")[0].split(".")
                product_smiles = rsmi.split(">")[-1]

                # Check if one reactant contains an epoxide
                epoxide_pattern = Chem.MolFromSmarts("[#6]1[#8][#6]1")
                phenol_pattern = Chem.MolFromSmarts("c[OH]")

                epoxide_found = False
                phenol_found = False

                for reactant in reactants_smiles:
                    try:
                        mol = Chem.MolFromSmiles(reactant)
                        if mol:
                            if mol.HasSubstructMatch(epoxide_pattern):
                                epoxide_found = True
                            if mol.HasSubstructMatch(phenol_pattern):
                                phenol_found = True
                    except:
                        continue

                # Check if product has a C-O-C linkage where one C is aromatic
                ether_pattern = Chem.MolFromSmarts("c-[#8]-[#6]")
                try:
                    product_mol = Chem.MolFromSmiles(product_smiles)
                    if (
                        product_mol
                        and product_mol.HasSubstructMatch(ether_pattern)
                        and epoxide_found
                        and phenol_found
                    ):
                        print("Found epoxide alkylation of phenol")
                        found_epoxide_alkylation = True
                except:
                    pass

        # Continue traversing
        for child in node.get("children", []):
            dfs_traverse(child)

    # Start traversal from the root
    dfs_traverse(route)

    return found_epoxide_alkylation
