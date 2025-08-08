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
    This function detects if stereochemistry is preserved throughout the synthesis.
    """
    all_steps_preserve_stereochemistry = True

    def dfs_traverse(node):
        nonlocal all_steps_preserve_stereochemistry

        if node["type"] == "reaction":
            # Extract reactants and product
            rsmi = node["metadata"]["rsmi"]
            reactants_smiles = rsmi.split(">")[0]
            product_smiles = rsmi.split(">")[-1]

            # Count stereocenters in reactants and products
            reactants_mol = Chem.MolFromSmiles(reactants_smiles)
            product_mol = Chem.MolFromSmiles(product_smiles)

            if reactants_mol and product_mol:
                # Get chiral centers
                reactant_chiral_centers = Chem.FindMolChiralCenters(
                    reactants_mol, includeUnassigned=True
                )
                product_chiral_centers = Chem.FindMolChiralCenters(
                    product_mol, includeUnassigned=True
                )

                # Simple check: if product has fewer stereocenters than reactants, stereochemistry might not be preserved
                # This is a simplified approach - would need refinement for production use
                if len(product_chiral_centers) < len(reactant_chiral_centers):
                    all_steps_preserve_stereochemistry = False
                    print("Potential loss of stereochemistry detected")

        # Process children
        for child in node.get("children", []):
            dfs_traverse(child)

    # Start traversal
    dfs_traverse(route)

    if all_steps_preserve_stereochemistry:
        print("Stereochemistry preservation confirmed")

    return all_steps_preserve_stereochemistry
