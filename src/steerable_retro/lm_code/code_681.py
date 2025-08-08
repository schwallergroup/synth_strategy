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
    This function detects if a synthetic route employs aromatic nitration,
    introducing a nitro group to an aromatic ring.
    """
    # Initialize tracking variable
    has_aromatic_nitration = False

    def dfs_traverse(node):
        nonlocal has_aromatic_nitration

        if node["type"] == "reaction":
            # Extract reactants and product from reaction SMILES
            rsmi = node.get("metadata", {}).get("rsmi", "")
            if not rsmi:
                return

            reactants_smiles = rsmi.split(">")[0].split(".")
            product_smiles = rsmi.split(">")[-1]

            try:
                reactants = [Chem.MolFromSmiles(r) for r in reactants_smiles if r]
                product = Chem.MolFromSmiles(product_smiles)

                if not all(reactants) or not product:
                    print("Failed to parse some molecules in reaction")
                    return

                # Check for aromatic ring in reactants without nitro
                aromatic_pattern = Chem.MolFromSmarts("c1ccccc1")
                nitro_pattern = Chem.MolFromSmarts("[#6]~[N+](=[O])[O-]")

                # Check if product has nitro group on aromatic ring
                if product.HasSubstructMatch(nitro_pattern):
                    # Check if any reactant is an aromatic compound
                    if any(r.HasSubstructMatch(aromatic_pattern) for r in reactants):
                        # Check if nitric acid or NO2+ source is among reactants
                        nitric_acid_pattern = Chem.MolFromSmarts("O[N+](=O)[O-]")
                        has_nitric_acid = any(
                            r.HasSubstructMatch(nitric_acid_pattern) for r in reactants if r
                        )

                        if has_nitric_acid:
                            print("Found aromatic nitration")
                            has_aromatic_nitration = True

            except Exception as e:
                print(f"Error processing reaction: {e}")

        # Recursively process children
        for child in node.get("children", []):
            dfs_traverse(child)

    # Start traversal from the root
    dfs_traverse(route)

    print(f"Aromatic nitration strategy detected: {has_aromatic_nitration}")
    return has_aromatic_nitration
