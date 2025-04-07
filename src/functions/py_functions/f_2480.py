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
    This function detects nucleophilic aromatic substitution where a chlorine
    on an aromatic ring is replaced by an amine.
    """
    nas_detected = False

    def dfs_traverse(node):
        nonlocal nas_detected

        if node["type"] == "reaction":
            if "metadata" in node and "rsmi" in node["metadata"]:
                rsmi = node["metadata"]["rsmi"]
                reactants = rsmi.split(">")[0].split(".")
                product = rsmi.split(">")[-1]

                # Check for nucleophilic aromatic substitution
                try:
                    # Pattern for chlorine on aromatic
                    cl_aromatic_pattern = Chem.MolFromSmarts("c-Cl")
                    # Pattern for amine on aromatic
                    amine_aromatic_pattern = Chem.MolFromSmarts("c-[N]")

                    # Check if any reactant has chlorine on aromatic
                    cl_in_reactants = False
                    for reactant in reactants:
                        react_mol = Chem.MolFromSmiles(reactant)
                        if react_mol and react_mol.HasSubstructMatch(
                            cl_aromatic_pattern
                        ):
                            cl_in_reactants = True
                            break

                    # Check if product has amine on aromatic where chlorine was
                    if cl_in_reactants:
                        prod_mol = Chem.MolFromSmiles(product)
                        if prod_mol and prod_mol.HasSubstructMatch(
                            amine_aromatic_pattern
                        ):
                            nas_detected = True
                            print("Nucleophilic aromatic substitution detected")
                except:
                    pass

        for child in node.get("children", []):
            dfs_traverse(child)

    dfs_traverse(route)
    return nas_detected
