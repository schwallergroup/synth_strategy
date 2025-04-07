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
    Detects an ester to acid to amide functional group transformation sequence.
    """
    # Track if we've found each step in the sequence
    found_ester_hydrolysis = False
    found_amide_formation = False

    def dfs_traverse(node):
        nonlocal found_ester_hydrolysis, found_amide_formation

        if (
            node["type"] == "reaction"
            and "metadata" in node
            and "rsmi" in node["metadata"]
        ):
            rsmi = node["metadata"]["rsmi"]
            reactants = rsmi.split(">")[0].split(".")
            product = rsmi.split(">")[-1]

            # Check for ester hydrolysis
            ester_hydrolysis_smarts = "[C:1](=[O:2])[O:3][C:4]>>[C:1](=[O:2])[O:3]"
            rxn_ester = AllChem.ReactionFromSmarts(ester_hydrolysis_smarts)

            # Check for amide formation
            amide_formation_smarts = "[C:1](=[O:2])[O,OH].[N:3]>>[C:1](=[O:2])[N:3]"
            rxn_amide = AllChem.ReactionFromSmarts(amide_formation_smarts)

            # Convert to RDKit molecules
            reactant_mols = [Chem.MolFromSmiles(r) for r in reactants]
            product_mol = Chem.MolFromSmiles(product)

            if product_mol and all(r for r in reactant_mols):
                # Check for ester hydrolysis
                if len(reactant_mols) == 1 and rxn_ester.RunReactants(
                    (reactant_mols[0],)
                ):
                    found_ester_hydrolysis = True
                    print("Found ester hydrolysis step")

                # Check for amide formation
                if len(reactant_mols) >= 2 and rxn_amide.RunReactants(
                    (reactant_mols[0], reactant_mols[1])
                ):
                    found_amide_formation = True
                    print("Found amide formation step")

        # Continue traversing
        for child in node.get("children", []):
            dfs_traverse(child)

    # Start traversal from the root
    dfs_traverse(route)

    # Return True if both steps were found
    return found_ester_hydrolysis and found_amide_formation
