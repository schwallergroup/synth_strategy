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
    This function detects the phenol → triflate → aryl-amine transformation sequence.
    """
    # Track transformations
    phenol_to_triflate = False
    triflate_to_arylamine = False

    # SMARTS patterns
    phenol_pattern = Chem.MolFromSmarts("c[OH]")
    triflate_pattern = Chem.MolFromSmarts("cOS(=O)(=O)C(F)(F)F")
    arylamine_pattern = Chem.MolFromSmarts("c[N]")

    def dfs_traverse(node):
        nonlocal phenol_to_triflate, triflate_to_arylamine

        if node["type"] == "reaction" and "rsmi" in node["metadata"]:
            rsmi = node["metadata"]["rsmi"]
            reactants_smiles = rsmi.split(">")[0]
            product_smiles = rsmi.split(">")[-1]

            reactants_mol = Chem.MolFromSmiles(reactants_smiles)
            product_mol = Chem.MolFromSmiles(product_smiles)

            if reactants_mol and product_mol:
                # Check for phenol → triflate transformation
                if (
                    reactants_mol.HasSubstructMatch(phenol_pattern)
                    and not reactants_mol.HasSubstructMatch(triflate_pattern)
                    and product_mol.HasSubstructMatch(triflate_pattern)
                ):
                    phenol_to_triflate = True
                    print("Detected phenol to triflate transformation")

                # Check for triflate → arylamine transformation
                if (
                    reactants_mol.HasSubstructMatch(triflate_pattern)
                    and product_mol.HasSubstructMatch(arylamine_pattern)
                    and not product_mol.HasSubstructMatch(triflate_pattern)
                ):
                    triflate_to_arylamine = True
                    print("Detected triflate to arylamine transformation")

        # Traverse children
        for child in node.get("children", []):
            dfs_traverse(child)

    # Start traversal
    dfs_traverse(route)

    # Strategy is present if we have both transformations
    strategy_present = phenol_to_triflate and triflate_to_arylamine

    if strategy_present:
        print("Phenol → triflate → aryl-amine transformation sequence detected")

    return strategy_present
