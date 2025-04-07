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
    This function detects late-stage coupling of a piperazine moiety.
    """
    late_stage_piperazine = False

    def dfs_traverse(node, depth=0):
        nonlocal late_stage_piperazine

        if (
            node["type"] == "reaction"
            and "metadata" in node
            and "rsmi" in node["metadata"]
        ):
            rsmi = node["metadata"]["rsmi"]
            reactants = rsmi.split(">")[0].split(".")
            product = rsmi.split(">")[-1]

            # Check if this is a late stage reaction (depth 0 or 1)
            if depth <= 1:
                # Check for piperazine in reactants
                piperazine_pattern = Chem.MolFromSmarts("[N]1CCN([C,H])CC1")

                for reactant in reactants:
                    try:
                        r_mol = Chem.MolFromSmiles(reactant)
                        if r_mol and r_mol.HasSubstructMatch(piperazine_pattern):
                            # Check if product has piperazine connected to aromatic
                            p_mol = Chem.MolFromSmiles(product)
                            if p_mol:
                                piperazine_aromatic_pattern = Chem.MolFromSmarts(
                                    "[n,c]~[N]1CCN([C,H])CC1"
                                )
                                if p_mol.HasSubstructMatch(piperazine_aromatic_pattern):
                                    print("Found late-stage piperazine coupling")
                                    late_stage_piperazine = True
                    except:
                        continue

        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    dfs_traverse(route)
    return late_stage_piperazine
