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
    This function detects a strategy involving the formation of a sulfonamide
    using a nitro-substituted sulfonyl chloride reagent.
    """
    # Initialize tracking variable
    has_nitro_sulfonamide_formation = False

    # SMARTS patterns
    nitro_pattern = Chem.MolFromSmarts("[#6]-[N+](=[O])[O-]")
    sulfonyl_chloride_pattern = Chem.MolFromSmarts("ClS(=O)(=O)[c]")

    def dfs_traverse(node):
        nonlocal has_nitro_sulfonamide_formation

        if node["type"] == "reaction" and "rsmi" in node.get("metadata", {}):
            rsmi = node["metadata"]["rsmi"]
            reactants_smiles = rsmi.split(">")[0]
            product_smiles = rsmi.split(">")[-1]

            reactants = [Chem.MolFromSmiles(r) for r in reactants_smiles.split(".")]
            product = Chem.MolFromSmiles(product_smiles)

            # Check for nitro-containing sulfonyl chloride
            nitro_sulfonyl_chloride_present = any(
                r
                and r.HasSubstructMatch(nitro_pattern)
                and r.HasSubstructMatch(sulfonyl_chloride_pattern)
                for r in reactants
                if r
            )

            # Check for sulfonamide formation
            if (
                nitro_sulfonyl_chloride_present
                and product
                and product.HasSubstructMatch(Chem.MolFromSmarts("[NH]S(=O)(=O)[c]"))
            ):
                print("Detected nitro-sulfonamide formation")
                has_nitro_sulfonamide_formation = True

        # Traverse children
        for child in node.get("children", []):
            dfs_traverse(child)

    # Start traversal
    dfs_traverse(route)

    return has_nitro_sulfonamide_formation
