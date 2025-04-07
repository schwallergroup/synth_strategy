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
    Detects a convergent synthesis where two complex fragments are combined via reductive amination.
    """
    reductive_amination_found = False
    late_stage = False

    def dfs_traverse(node, depth=0):
        nonlocal reductive_amination_found, late_stage

        if node["type"] == "reaction":
            # Check if this is a reductive amination reaction
            if "rsmi" in node.get("metadata", {}):
                rsmi = node["metadata"]["rsmi"]
                reactants = rsmi.split(">")[0].split(".")
                product = rsmi.split(">")[-1]

                # Check if we have at least 2 reactants (convergent)
                if len(reactants) >= 2:
                    # Check for aldehyde in one reactant
                    aldehyde_pattern = Chem.MolFromSmarts("[CH]=[O]")
                    # Check for amine in another reactant
                    amine_pattern = Chem.MolFromSmarts("[NH2]")
                    # Check for secondary amine in product
                    sec_amine_pattern = Chem.MolFromSmarts("[NH]")

                    aldehyde_found = False
                    amine_found = False

                    for reactant in reactants:
                        mol = Chem.MolFromSmiles(reactant)
                        if mol:
                            if mol.HasSubstructMatch(aldehyde_pattern):
                                aldehyde_found = True
                            if mol.HasSubstructMatch(amine_pattern):
                                amine_found = True

                    product_mol = Chem.MolFromSmiles(product)
                    sec_amine_found = False
                    if product_mol and product_mol.HasSubstructMatch(sec_amine_pattern):
                        sec_amine_found = True

                    if aldehyde_found and amine_found and sec_amine_found:
                        reductive_amination_found = True
                        # If depth is 0 or 1, it's considered late-stage
                        if depth <= 1:
                            late_stage = True
                            print(f"Found late-stage reductive amination at depth {depth}")

        # Traverse children
        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    dfs_traverse(route)

    if reductive_amination_found and late_stage:
        print("Detected convergent synthesis with late-stage reductive amination")
        return True
    return False
