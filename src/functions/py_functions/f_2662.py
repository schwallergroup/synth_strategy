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
    This function detects a specific cascade of functional group transformations:
    ketone → alcohol → ether → amine → amide
    """
    # Initialize tracking variables
    transformations = {
        "ketone_to_alcohol": False,
        "alcohol_to_ether": False,
        "nitrile_to_amine": False,
        "amine_to_amide": False,
    }

    def dfs_traverse(node):
        if node["type"] == "reaction":
            if "rsmi" in node.get("metadata", {}):
                rsmi = node["metadata"]["rsmi"]
                reactants = rsmi.split(">")[0].split(".")
                product = rsmi.split(">")[-1]

                # Convert to RDKit molecules
                reactant_mols = [
                    Chem.MolFromSmiles(r) for r in reactants if Chem.MolFromSmiles(r)
                ]
                product_mol = Chem.MolFromSmiles(product)

                if product_mol and reactant_mols:
                    # Check for ketone to alcohol reduction
                    ketone_pattern = Chem.MolFromSmarts("[C](=[O])[#6]")
                    alcohol_pattern = Chem.MolFromSmarts("[C]([OH])[#6]")

                    if any(
                        r.HasSubstructMatch(ketone_pattern) for r in reactant_mols if r
                    ) and product_mol.HasSubstructMatch(alcohol_pattern):
                        print("Found ketone to alcohol transformation")
                        transformations["ketone_to_alcohol"] = True

                    # Check for alcohol to ether conversion
                    if any(
                        r.HasSubstructMatch(alcohol_pattern) for r in reactant_mols if r
                    ):
                        ether_pattern = Chem.MolFromSmarts("[C]([O][C])[#6]")
                        if product_mol.HasSubstructMatch(ether_pattern):
                            print("Found alcohol to ether transformation")
                            transformations["alcohol_to_ether"] = True

                    # Check for nitrile to amine reduction
                    nitrile_pattern = Chem.MolFromSmarts("[C]#[N]")
                    amine_pattern = Chem.MolFromSmarts("[C][NH2]")

                    if any(
                        r.HasSubstructMatch(nitrile_pattern) for r in reactant_mols if r
                    ) and product_mol.HasSubstructMatch(amine_pattern):
                        print("Found nitrile to amine transformation")
                        transformations["nitrile_to_amine"] = True

                    # Check for amine to amide formation
                    if any(
                        r.HasSubstructMatch(amine_pattern) for r in reactant_mols if r
                    ):
                        amide_pattern = Chem.MolFromSmarts("[C](=[O])[N]")
                        if product_mol.HasSubstructMatch(amide_pattern):
                            print("Found amine to amide transformation")
                            transformations["amine_to_amide"] = True

        # Traverse children
        for child in node.get("children", []):
            dfs_traverse(child)

    # Start traversal
    dfs_traverse(route)

    # Check if we have at least 3 of the 4 transformations
    transformation_count = sum(transformations.values())
    return transformation_count >= 3
