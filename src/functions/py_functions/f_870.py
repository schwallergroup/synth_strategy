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
    This function detects if the synthesis includes a nitro reduction step
    (converting -NO2 to -NH2 on an aromatic ring).
    """
    nitro_reduction_found = False

    def dfs_traverse(node):
        nonlocal nitro_reduction_found

        if node["type"] == "reaction" and node.get("metadata", {}).get("rsmi"):
            rsmi = node["metadata"]["rsmi"]
            reactants_part = rsmi.split(">")[0]
            product_part = rsmi.split(">")[-1]

            # Pattern for nitro group
            nitro_pattern = Chem.MolFromSmarts("[c][N+](=[O])[O-]")
            # Pattern for amine group
            amine_pattern = Chem.MolFromSmarts("[c][NH2]")

            try:
                reactant_mols = [
                    Chem.MolFromSmiles(r) for r in reactants_part.split(".")
                ]
                product_mol = Chem.MolFromSmiles(product_part)

                # Check for nitro reduction
                if (
                    any(
                        r and r.HasSubstructMatch(nitro_pattern)
                        for r in reactant_mols
                        if r
                    )
                    and product_mol
                    and product_mol.HasSubstructMatch(amine_pattern)
                ):
                    nitro_reduction_found = True
                    print("Detected nitro reduction step")
            except:
                pass

        # Continue traversing
        for child in node.get("children", []):
            dfs_traverse(child)

    # Start traversal from the root
    dfs_traverse(route)
    return nitro_reduction_found
