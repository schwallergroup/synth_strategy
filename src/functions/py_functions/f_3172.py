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
    This function detects if the synthetic route involves Weinreb amide chemistry
    (formation and subsequent transformation).
    """
    weinreb_formation = False
    weinreb_transformation = False

    # SMARTS pattern for Weinreb amide
    weinreb_pattern = Chem.MolFromSmarts("[#6]-[#7](-[#8][#6])-[#6]=[#8]")

    def dfs_traverse(node):
        nonlocal weinreb_formation, weinreb_transformation

        if (
            node["type"] == "reaction"
            and "metadata" in node
            and "rsmi" in node["metadata"]
        ):
            rsmi = node["metadata"]["rsmi"]
            reactants = rsmi.split(">")[0].split(".")
            product = rsmi.split(">")[-1]

            # Check for Weinreb amide formation
            product_mol = Chem.MolFromSmiles(product)
            if product_mol and product_mol.HasSubstructMatch(weinreb_pattern):
                # Check if Weinreb amide was not in reactants
                weinreb_in_reactants = False
                for reactant in reactants:
                    reactant_mol = Chem.MolFromSmiles(reactant)
                    if reactant_mol and reactant_mol.HasSubstructMatch(weinreb_pattern):
                        weinreb_in_reactants = True
                        break

                if not weinreb_in_reactants:
                    print("Weinreb amide formation detected in reaction:", rsmi)
                    weinreb_formation = True

            # Check for Weinreb amide transformation
            weinreb_in_reactants = False
            for reactant in reactants:
                reactant_mol = Chem.MolFromSmiles(reactant)
                if reactant_mol and reactant_mol.HasSubstructMatch(weinreb_pattern):
                    weinreb_in_reactants = True
                    break

            if weinreb_in_reactants:
                product_mol = Chem.MolFromSmiles(product)
                if product_mol and not product_mol.HasSubstructMatch(weinreb_pattern):
                    print("Weinreb amide transformation detected in reaction:", rsmi)
                    weinreb_transformation = True

        # Continue traversing
        for child in node.get("children", []):
            dfs_traverse(child)

    # Start traversal
    dfs_traverse(route)
    return weinreb_formation and weinreb_transformation
