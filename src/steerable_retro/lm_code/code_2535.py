#!/bin/python

"""LM-defined function for strategy description."""

from rdkit.Chem import AllChem, rdFMCS
import copy
from collections import deque
import rdkit.Chem as Chem
from rdkit.Chem import rdMolDescriptors
from rdkit.Chem import rdChemReactions
from rdkit.Chem import AllChem
from rdkit.Chem import rdFMCS
import rdkit.Chem.rdFMCS
from rdkit.Chem import AllChem, Descriptors, rdMolDescriptors
from rdkit import Chem
from rdkit.Chem import Descriptors
from rdkit.Chem import AllChem, rdMolDescriptors
from rdkit.Chem import AllChem, Descriptors, Lipinski
from rdkit.Chem import rdmolops
import re
from rdkit.Chem.Scaffolds import MurckoScaffold
from rdkit.Chem import AllChem, Descriptors
import traceback
import rdkit
from collections import Counter


def main(route):
    """
    Detects a strategy involving benzylic functionalization, specifically
    the transformation of a methyl group to a bromomethyl group, followed by
    N-alkylation to form a benzylic amine.
    """
    # Track if we've seen each transformation
    methyl_to_bromomethyl = False
    bromomethyl_to_benzylamine = False

    def dfs_traverse(node):
        nonlocal methyl_to_bromomethyl, bromomethyl_to_benzylamine

        if node["type"] == "reaction":
            rsmi = node["metadata"].get("rsmi", "")

            # Extract reactants and product
            parts = rsmi.split(">")
            if len(parts) >= 3:
                reactants = parts[0].split(".")
                product = parts[2]

                # Check for methyl to bromomethyl conversion
                methyl_pattern = Chem.MolFromSmarts("[c]-[CH3]")
                bromomethyl_pattern = Chem.MolFromSmarts("[c]-[CH2][Br]")

                for reactant in reactants:
                    reactant_mol = Chem.MolFromSmiles(reactant)
                    if reactant_mol and reactant_mol.HasSubstructMatch(methyl_pattern):
                        product_mol = Chem.MolFromSmiles(product)
                        if product_mol and product_mol.HasSubstructMatch(bromomethyl_pattern):
                            print("Detected methyl to bromomethyl conversion")
                            methyl_to_bromomethyl = True

                # Check for bromomethyl to benzylamine conversion
                bromomethyl_pattern = Chem.MolFromSmarts("[c]-[CH2][Br]")
                benzylamine_pattern = Chem.MolFromSmarts("[c]-[CH2]-[N]")

                for reactant in reactants:
                    reactant_mol = Chem.MolFromSmiles(reactant)
                    if reactant_mol and reactant_mol.HasSubstructMatch(bromomethyl_pattern):
                        product_mol = Chem.MolFromSmiles(product)
                        if product_mol and product_mol.HasSubstructMatch(benzylamine_pattern):
                            print("Detected bromomethyl to benzylamine conversion")
                            bromomethyl_to_benzylamine = True

        # Continue traversing
        for child in node.get("children", []):
            dfs_traverse(child)

    # Start traversal from the root
    dfs_traverse(route)

    # Return True if both transformations are found
    return methyl_to_bromomethyl and bromomethyl_to_benzylamine
