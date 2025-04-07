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
    Detects if the synthesis route involves installation of a difluoromethoxy group.
    """
    difluoromethoxy_detected = False

    def dfs_traverse(node, depth=0):
        nonlocal difluoromethoxy_detected

        if node["type"] == "reaction":
            rsmi = node["metadata"]["rsmi"]
            reactants_smiles = rsmi.split(">")[0].split(".")
            product_smiles = rsmi.split(">")[-1]

            # Create RDKit molecules
            try:
                reactant_mols = [Chem.MolFromSmiles(r) for r in reactants_smiles]
                product_mol = Chem.MolFromSmiles(product_smiles)

                if product_mol and all(r for r in reactant_mols):
                    # Check for difluoromethoxy group in product
                    difluoromethoxy_pattern = Chem.MolFromSmarts("[O][C]([F])[F]")

                    # Check if pattern exists in product but not in reactants
                    if product_mol.HasSubstructMatch(difluoromethoxy_pattern):
                        difluoromethoxy_in_reactants = False
                        for r_mol in reactant_mols:
                            if r_mol.HasSubstructMatch(difluoromethoxy_pattern):
                                difluoromethoxy_in_reactants = True
                                break

                        if not difluoromethoxy_in_reactants:
                            difluoromethoxy_detected = True
                            print(
                                f"Difluoromethoxy installation detected at depth {depth}"
                            )
                            return
            except:
                print(f"Error processing reaction at depth {depth}")

        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    dfs_traverse(route)
    return difluoromethoxy_detected
