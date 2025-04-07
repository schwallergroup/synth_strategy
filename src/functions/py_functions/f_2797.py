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
    This function detects if the synthetic route involves a biaryl coupling strategy.
    It looks for a C-C bond formation between two aromatic rings.
    """
    biaryl_coupling_found = False

    def dfs_traverse(node):
        nonlocal biaryl_coupling_found

        if (
            node["type"] == "reaction"
            and "metadata" in node
            and "rsmi" in node["metadata"]
        ):
            rsmi = node["metadata"]["rsmi"]
            reactants_part = rsmi.split(">")[0]
            products_part = rsmi.split(">")[-1]

            # Check if this is a biaryl coupling reaction
            reactants = reactants_part.split(".")

            # Look for reactions where we have aromatic rings in the reactants
            # and a new biaryl bond in the product
            if len(reactants) >= 2:
                try:
                    product_mol = Chem.MolFromSmiles(products_part)
                    reactant_mols = [Chem.MolFromSmiles(r) for r in reactants]

                    # Check if product has a biaryl bond that's not in any reactant
                    if product_mol and all(reactant_mols):
                        # Simple heuristic: if we have two aromatic reactants and
                        # the product contains both aromatic systems connected
                        aromatic_reactants = 0
                        for r_mol in reactant_mols:
                            if any(atom.GetIsAromatic() for atom in r_mol.GetAtoms()):
                                aromatic_reactants += 1

                        if aromatic_reactants >= 2 and any(
                            atom.GetIsAromatic() for atom in product_mol.GetAtoms()
                        ):
                            print("Potential biaryl coupling detected")
                            biaryl_coupling_found = True
                except:
                    pass

        for child in node.get("children", []):
            dfs_traverse(child)

    dfs_traverse(route)
    return biaryl_coupling_found
