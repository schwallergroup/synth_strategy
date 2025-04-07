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
    Detects if the synthesis route includes phenol O-methylation.
    Looks for conversion of phenolic OH to methyl ether.
    """
    methylation_count = 0

    def dfs_traverse(node):
        nonlocal methylation_count

        if node["type"] == "reaction" and "metadata" in node and "rsmi" in node["metadata"]:
            rsmi = node["metadata"]["rsmi"]
            reactants = rsmi.split(">")[0].split(".")
            product = rsmi.split(">")[-1]

            # Check for phenol to methyl ether conversion
            for reactant in reactants:
                reactant_mol = Chem.MolFromSmiles(reactant)
                product_mol = Chem.MolFromSmiles(product)

                if reactant_mol and product_mol:
                    phenol_patt = Chem.MolFromSmarts("[c][OH]")
                    methoxy_patt = Chem.MolFromSmarts("[c][O][CH3]")

                    if (
                        reactant_mol.HasSubstructMatch(phenol_patt)
                        and product_mol.HasSubstructMatch(methoxy_patt)
                        and not reactant_mol.HasSubstructMatch(methoxy_patt)
                    ):
                        print(f"Phenol methylation detected: {rsmi}")
                        methylation_count += 1

        for child in node.get("children", []):
            dfs_traverse(child)

    dfs_traverse(route)
    return methylation_count >= 1
