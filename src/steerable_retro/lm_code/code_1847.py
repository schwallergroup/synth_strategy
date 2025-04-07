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
    This function detects TBDMS protection of an alcohol in the synthetic route.
    """
    tbdms_protection_found = False

    def dfs_traverse(node):
        nonlocal tbdms_protection_found

        if node["type"] == "reaction":
            if "rsmi" in node.get("metadata", {}):
                rsmi = node["metadata"]["rsmi"]
                reactants = rsmi.split(">")[0].split(".")
                product = rsmi.split(">")[-1]

                # Check for TBDMS protection
                alcohol_pattern = Chem.MolFromSmarts("[#6]-[#8;H1]")
                tbdms_cl_pattern = Chem.MolFromSmarts("Cl[Si](C)(C)C(C)(C)C")
                protected_alcohol_pattern = Chem.MolFromSmarts("[#6]-[#8]-[Si](C)(C)C(C)(C)C")

                try:
                    product_mol = Chem.MolFromSmiles(product)
                    if (
                        product_mol
                        and alcohol_pattern
                        and tbdms_cl_pattern
                        and protected_alcohol_pattern
                    ):
                        for reactant in reactants:
                            reactant_mol = Chem.MolFromSmiles(reactant)
                            if reactant_mol and reactant_mol.HasSubstructMatch(alcohol_pattern):
                                for r in reactants:
                                    r_mol = Chem.MolFromSmiles(r)
                                    if r_mol and r_mol.HasSubstructMatch(tbdms_cl_pattern):
                                        if product_mol.HasSubstructMatch(protected_alcohol_pattern):
                                            print("Found TBDMS protection of alcohol")
                                            tbdms_protection_found = True
                except:
                    print("Error in SMILES processing for TBDMS protection detection")

        for child in node.get("children", []):
            dfs_traverse(child)

    dfs_traverse(route)
    return tbdms_protection_found
