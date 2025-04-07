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
    This function detects if the synthesis uses a convergent strategy with
    amide coupling to join two complex fragments.
    """
    has_convergent_amide_coupling = False

    def dfs_traverse(node, depth=0):
        nonlocal has_convergent_amide_coupling

        if node["type"] == "reaction":
            if "rsmi" in node["metadata"]:
                rsmi = node["metadata"]["rsmi"]
                reactants = rsmi.split(">")[0].split(".")
                product = rsmi.split(">")[-1]

                # Check if this is an amide formation
                product_mol = Chem.MolFromSmiles(product)
                amide_pattern = Chem.MolFromSmarts("[#7]-[#6](=[#8])-[#6]")

                # Check if this is a convergent step (has multiple complex reactants)
                if (
                    product_mol
                    and product_mol.HasSubstructMatch(amide_pattern)
                    and len(reactants) >= 2
                ):
                    # Check complexity of reactants (both should have aromatic rings)
                    complex_reactants = 0
                    for reactant in reactants:
                        reactant_mol = Chem.MolFromSmiles(reactant)
                        if reactant_mol:
                            aromatic_pattern = Chem.MolFromSmarts("a")
                            if reactant_mol.HasSubstructMatch(aromatic_pattern):
                                complex_reactants += 1

                    if complex_reactants >= 2:
                        print(f"Found convergent amide coupling at depth {depth}")
                        has_convergent_amide_coupling = True

        # Traverse children
        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    # Start traversal from root
    dfs_traverse(route)
    return has_convergent_amide_coupling
