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
    This function detects a strategy involving nitro group removal.
    """
    # Track if we found nitro removal
    found_nitro_removal = False

    # Define SMARTS pattern for nitro group
    nitro_pattern = Chem.MolFromSmarts("[c][N+](=[O])[O-]")

    def dfs_traverse(node, depth=0):
        nonlocal found_nitro_removal

        if node["type"] == "reaction":
            # Extract reactants and product
            rsmi = node["metadata"]["rsmi"]
            reactants_smiles = rsmi.split(">")[0].split(".")
            product_smiles = rsmi.split(">")[-1].split(".")

            try:
                # Check if reactant has nitro but product doesn't
                reactant_mols = [Chem.MolFromSmiles(r) for r in reactants_smiles]
                product_mols = [Chem.MolFromSmiles(p) for p in product_smiles]

                reactant_has_nitro = any(
                    r and r.HasSubstructMatch(nitro_pattern) for r in reactant_mols
                )

                # Check if main product doesn't have nitro
                # We assume the first product is the main one
                if product_mols and len(product_mols) > 0 and product_mols[0]:
                    product_has_nitro = product_mols[0].HasSubstructMatch(nitro_pattern)

                    if reactant_has_nitro and not product_has_nitro:
                        found_nitro_removal = True
                        print(f"Found nitro removal at depth {depth}")

            except:
                print(f"Error processing SMILES at depth {depth}")

        # Process children
        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    # Start traversal
    dfs_traverse(route)

    print(f"Nitro removal strategy detected: {found_nitro_removal}")
    return found_nitro_removal
