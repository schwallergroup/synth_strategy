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
    Detects if the synthesis involves sequential functionalization of aromatic/heteroaromatic systems.
    """
    functionalization_steps = 0

    def dfs_traverse(node):
        nonlocal functionalization_steps

        if node["type"] == "reaction":
            if "rsmi" in node.get("metadata", {}):
                rsmi = node["metadata"]["rsmi"]
                product = rsmi.split(">")[-1]
                reactants = rsmi.split(">")[0].split(".")

                # Check if this is an aromatic functionalization
                product_mol = Chem.MolFromSmiles(product)
                if (
                    product_mol and product_mol.GetNumAtoms() > 5
                ):  # Ensure it's not a trivial molecule
                    aromatic_pattern = Chem.MolFromSmarts("c")
                    if product_mol.HasSubstructMatch(aromatic_pattern):
                        # Common functionalization patterns
                        patterns = [
                            ("Halogenation", Chem.MolFromSmarts("c[Br,I,Cl]")),
                            ("Amination", Chem.MolFromSmarts("c[N]")),
                            ("Alkoxylation", Chem.MolFromSmarts("c[O][C]")),
                            ("Metallation", Chem.MolFromSmarts("c[Sn,B,Zn,Mg]")),
                            ("Sulfonylation", Chem.MolFromSmarts("c[S](=[O])(=[O])[O,N]")),
                        ]

                        for name, pattern in patterns:
                            if product_mol.HasSubstructMatch(pattern):
                                # Check if this is a new functionality
                                is_new = True
                                for reactant in reactants:
                                    reactant_mol = Chem.MolFromSmiles(reactant)
                                    if reactant_mol and reactant_mol.HasSubstructMatch(pattern):
                                        is_new = False
                                        break

                                if is_new:
                                    print(f"Found {name} of aromatic system")
                                    functionalization_steps += 1
                                    break

        # Traverse children
        for child in node.get("children", []):
            dfs_traverse(child)

    dfs_traverse(route)
    return functionalization_steps >= 2  # At least 2 functionalization steps
