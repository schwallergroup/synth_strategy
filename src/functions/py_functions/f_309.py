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
    This function detects if stereochemistry is preserved throughout the synthesis.
    """
    has_stereocenter = False
    stereocenter_preserved = True

    def dfs_traverse(node):
        nonlocal has_stereocenter, stereocenter_preserved

        if node["type"] == "mol" and "smiles" in node:
            try:
                mol = Chem.MolFromSmiles(node["smiles"])
                if mol:
                    # Check for presence of stereochemistry in the molecule
                    chiral_centers = Chem.FindMolChiralCenters(
                        mol, includeUnassigned=False
                    )
                    if chiral_centers:
                        has_stereocenter = True
                        print(f"Detected stereocenter in molecule: {node['smiles']}")
            except:
                pass

        elif (
            node["type"] == "reaction"
            and "metadata" in node
            and "rsmi" in node["metadata"]
        ):
            rsmi = node["metadata"]["rsmi"]
            reactants_part = rsmi.split(">")[0]
            product_part = rsmi.split(">")[-1]

            try:
                # Check if stereochemistry is preserved in this reaction
                reactants = reactants_part.split(".")
                product_mol = Chem.MolFromSmiles(product_part)

                # Check for stereocenters in reactants
                reactant_has_stereocenter = False
                for reactant in reactants:
                    reactant_mol = Chem.MolFromSmiles(reactant)
                    if reactant_mol:
                        chiral_centers = Chem.FindMolChiralCenters(
                            reactant_mol, includeUnassigned=False
                        )
                        if chiral_centers:
                            reactant_has_stereocenter = True
                            break

                # Check for stereocenters in product
                product_has_stereocenter = False
                if product_mol:
                    chiral_centers = Chem.FindMolChiralCenters(
                        product_mol, includeUnassigned=False
                    )
                    if chiral_centers:
                        product_has_stereocenter = True

                # If reactant had stereocenter but product doesn't, stereochemistry was lost
                if reactant_has_stereocenter and not product_has_stereocenter:
                    stereocenter_preserved = False
                    print(f"Stereochemistry lost in reaction: {rsmi}")
            except:
                pass

        # Continue traversing
        for child in node.get("children", []):
            dfs_traverse(child)

    # Start traversal from the root
    dfs_traverse(route)

    # Return True if we found stereocenters and they were preserved
    return has_stereocenter and stereocenter_preserved
