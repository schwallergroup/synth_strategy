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
    Detects the use of Williamson ether synthesis for aryl methyl ether formation
    as part of the overall synthetic strategy.
    """
    found_ether_synthesis = False

    def dfs_traverse(node):
        nonlocal found_ether_synthesis

        if node["type"] == "reaction":
            if "metadata" in node and "rsmi" in node["metadata"]:
                rsmi = node["metadata"]["rsmi"]
                reactants = rsmi.split(">")[0].split(".")
                product = rsmi.split(">")[-1]

                # Check for Williamson ether synthesis
                if len(reactants) == 2:
                    phenol_pattern = Chem.MolFromSmarts("[OH][c]")
                    methyl_halide_pattern = Chem.MolFromSmarts("[CH3][I,Br,Cl,F]")
                    methoxy_pattern = Chem.MolFromSmarts("[CH3][O][c]")

                    reactant_mols = [Chem.MolFromSmiles(r) for r in reactants]
                    product_mol = Chem.MolFromSmiles(product)

                    if (
                        product_mol
                        and any(
                            m and m.HasSubstructMatch(phenol_pattern)
                            for m in reactant_mols
                        )
                        and any(
                            m and m.HasSubstructMatch(methyl_halide_pattern)
                            for m in reactant_mols
                        )
                        and product_mol.HasSubstructMatch(methoxy_pattern)
                    ):
                        print("Found Williamson ether synthesis")
                        found_ether_synthesis = True

        # Traverse children
        for child in node.get("children", []):
            dfs_traverse(child)

    # Start traversal
    dfs_traverse(route)

    return found_ether_synthesis
