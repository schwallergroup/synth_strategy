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
    Detects biaryl formation via Suzuki coupling.
    Looks for C-C bond formation between two aromatic rings where one reactant
    contains a boronic acid/ester and the other contains a halogen.
    """
    has_suzuki_coupling = False

    def dfs_traverse(node, depth=0):
        nonlocal has_suzuki_coupling

        if node["type"] == "reaction":
            if "rsmi" in node["metadata"]:
                rsmi = node["metadata"]["rsmi"]
                reactants = rsmi.split(">")[0].split(".")
                product = rsmi.split(">")[-1]

                # Check for Suzuki coupling pattern
                if len(reactants) >= 2:
                    try:
                        product_mol = Chem.MolFromSmiles(product)

                        # Check if product has biaryl system
                        biaryl_pattern = Chem.MolFromSmarts("c-c")
                        if product_mol.HasSubstructMatch(biaryl_pattern):

                            # Check if one reactant has boronic acid/ester
                            boronic_pattern = Chem.MolFromSmarts(
                                "[#6]-[#5](-[#8])-[#8]"
                            )

                            # Check if another reactant has halogen on aromatic
                            halo_aromatic_pattern = Chem.MolFromSmarts(
                                "c-[#9,#17,#35,#53]"
                            )

                            has_boronic = False
                            has_halo = False

                            for r in reactants:
                                r_mol = Chem.MolFromSmiles(r)
                                if r_mol:
                                    if r_mol.HasSubstructMatch(boronic_pattern):
                                        has_boronic = True
                                    if r_mol.HasSubstructMatch(halo_aromatic_pattern):
                                        has_halo = True

                            if has_boronic and has_halo:
                                has_suzuki_coupling = True
                                print(f"Found Suzuki coupling at depth {depth}")
                    except:
                        pass

        # Traverse children
        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    dfs_traverse(route)
    return has_suzuki_coupling
