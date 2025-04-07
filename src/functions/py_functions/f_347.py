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
    This function detects a linear synthesis strategy that maintains stereocenters
    throughout the synthesis.
    """
    linear_synthesis = True
    has_stereocenters = False

    def dfs_traverse(node):
        nonlocal linear_synthesis, has_stereocenters

        if node["type"] == "mol" and "smiles" in node:
            smiles = node["smiles"]

            try:
                mol = Chem.MolFromSmiles(smiles)
                if mol:
                    # Check for stereocenters
                    chiral_centers = Chem.FindMolChiralCenters(mol)
                    if len(chiral_centers) > 0:
                        has_stereocenters = True
                        print(
                            f"Detected {len(chiral_centers)} stereocenters in molecule: {smiles}"
                        )
            except:
                pass

        # Check if this is a convergent synthesis (more than 2 reactants)
        if (
            node["type"] == "reaction"
            and "metadata" in node
            and "rsmi" in node["metadata"]
        ):
            rsmi = node["metadata"]["rsmi"]
            reactants = rsmi.split(">")[0].split(".")
            if len(reactants) > 2:
                linear_synthesis = False
                print("Detected convergent step with more than 2 reactants:", rsmi)

        children = node.get("children", [])
        if len(children) > 2:
            linear_synthesis = False
            print("Detected convergent synthesis pattern with multiple children")

        for child in children:
            dfs_traverse(child)

    dfs_traverse(route)

    if linear_synthesis and has_stereocenters:
        print("Detected linear synthesis strategy that maintains stereocenters")

    return linear_synthesis and has_stereocenters
