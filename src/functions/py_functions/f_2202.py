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
    This function detects if the synthesis route involves biaryl coupling.
    """
    biaryl_coupling_found = False

    def dfs_traverse(node):
        nonlocal biaryl_coupling_found

        if node["type"] == "reaction" and "rsmi" in node.get("metadata", {}):
            rsmi = node["metadata"]["rsmi"]
            reactants = rsmi.split(">")[0].split(".")
            product = rsmi.split(">")[-1]

            # Check if product contains biaryl but reactants don't
            product_mol = Chem.MolFromSmiles(product)
            if product_mol:
                biaryl_pattern = Chem.MolFromSmarts("[c]-[c]")
                if product_mol.HasSubstructMatch(biaryl_pattern):
                    # Check if any single reactant has the biaryl pattern
                    biaryl_in_single_reactant = False
                    for reactant in reactants:
                        reactant_mol = Chem.MolFromSmiles(reactant)
                        if reactant_mol and reactant_mol.HasSubstructMatch(
                            biaryl_pattern
                        ):
                            biaryl_in_single_reactant = True

                    # If biaryl is in product but not in any single reactant, it's a biaryl coupling
                    if not biaryl_in_single_reactant and len(reactants) >= 2:
                        # Check if one reactant has a triflate group (common for biaryl coupling)
                        triflate_in_reactants = False
                        for reactant in reactants:
                            reactant_mol = Chem.MolFromSmiles(reactant)
                            if reactant_mol:
                                triflate_pattern = Chem.MolFromSmarts(
                                    "[#8]-[#16](=[#8])(=[#8])-[#6]([F])([F])[F]"
                                )
                                if reactant_mol.HasSubstructMatch(triflate_pattern):
                                    triflate_in_reactants = True

                        if triflate_in_reactants:
                            biaryl_coupling_found = True
                            print("Biaryl coupling with triflate detected")

        for child in node.get("children", []):
            dfs_traverse(child)

    dfs_traverse(route)
    return biaryl_coupling_found
