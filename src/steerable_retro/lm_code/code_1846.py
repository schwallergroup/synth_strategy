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
    This function detects N-alkylation of indole with a dibromoalkane followed by
    functionalization of the alkyl chain.
    """
    n_alkylation_found = False

    def dfs_traverse(node):
        nonlocal n_alkylation_found

        if node["type"] == "reaction":
            if "rsmi" in node.get("metadata", {}):
                rsmi = node["metadata"]["rsmi"]
                reactants = rsmi.split(">")[0].split(".")
                product = rsmi.split(">")[-1]

                # Check for indole N-alkylation with dibromoalkane
                indole_pattern = Chem.MolFromSmarts("c1ccc2[nH]ccc2c1")
                dibromo_pattern = Chem.MolFromSmarts("BrCCCBr")
                n_alkylated_indole_pattern = Chem.MolFromSmarts("c1ccc2n(CCC*)ccc2c1")

                try:
                    product_mol = Chem.MolFromSmiles(product)
                    if (
                        product_mol
                        and indole_pattern
                        and dibromo_pattern
                        and n_alkylated_indole_pattern
                    ):
                        for reactant in reactants:
                            reactant_mol = Chem.MolFromSmiles(reactant)
                            if reactant_mol and reactant_mol.HasSubstructMatch(indole_pattern):
                                for r in reactants:
                                    r_mol = Chem.MolFromSmiles(r)
                                    if r_mol and r_mol.HasSubstructMatch(dibromo_pattern):
                                        if product_mol.HasSubstructMatch(
                                            n_alkylated_indole_pattern
                                        ):
                                            print("Found indole N-alkylation with dibromoalkane")
                                            n_alkylation_found = True
                except:
                    print("Error in SMILES processing for N-alkylation detection")

        for child in node.get("children", []):
            dfs_traverse(child)

    dfs_traverse(route)
    return n_alkylation_found
