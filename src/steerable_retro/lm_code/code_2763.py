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
    This function detects if the synthetic route preserves stereocenters throughout the synthesis.
    """
    has_stereocenter = False
    preserves_stereocenter = False

    def dfs_traverse(node):
        nonlocal has_stereocenter, preserves_stereocenter

        if node["type"] == "mol":
            # Check if the final molecule has a stereocenter
            mol = Chem.MolFromSmiles(node["smiles"])
            if mol and Chem.FindMolChiralCenters(mol):
                has_stereocenter = True
                print(f"Found stereocenter in molecule: {node['smiles']}")

        elif node["type"] == "reaction":
            # Extract reactants and product
            rsmi = node["metadata"]["rsmi"]
            reactants_smiles = rsmi.split(">")[0].split(".")
            product_smiles = rsmi.split(">")[-1]

            product_mol = Chem.MolFromSmiles(product_smiles)
            product_chiral_centers = Chem.FindMolChiralCenters(product_mol) if product_mol else []

            if product_chiral_centers:
                # Check if any reactant also has the same chiral center
                for reactant in reactants_smiles:
                    reactant_mol = Chem.MolFromSmiles(reactant)
                    if reactant_mol and Chem.FindMolChiralCenters(reactant_mol):
                        preserves_stereocenter = True
                        print(
                            f"Found stereocenter preservation at depth {node.get('depth', 'unknown')}"
                        )
                        break

        # Traverse children
        for child in node.get("children", []):
            dfs_traverse(child)

    # Start traversal
    dfs_traverse(route)

    # Return True if the route has stereocenters and preserves them
    result = has_stereocenter and preserves_stereocenter
    print(f"Stereocenter preservation strategy detected: {result}")
    return result
