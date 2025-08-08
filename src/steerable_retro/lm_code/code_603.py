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
    This function detects if the synthetic route involves benzyl ether formation
    in a linear synthesis strategy.
    """
    # Track if we found benzyl ether formation
    benzyl_ether_formation_found = False

    def dfs_traverse(node):
        nonlocal benzyl_ether_formation_found

        # Only process reaction nodes
        if node["type"] == "reaction" and "metadata" in node and "rsmi" in node["metadata"]:
            rsmi = node["metadata"]["rsmi"]
            reactants_smiles = rsmi.split(">")[0].split(".")
            product_smiles = rsmi.split(">")[-1]

            # Convert SMILES to RDKit molecules
            reactants = [Chem.MolFromSmiles(r) for r in reactants_smiles if r]
            product = Chem.MolFromSmiles(product_smiles) if product_smiles else None

            if product and all(r is not None for r in reactants):
                # Check if any reactant contains phenol
                phenol_pattern = Chem.MolFromSmarts("c[OH]")
                has_phenol = any(r.HasSubstructMatch(phenol_pattern) for r in reactants)

                # Check if any reactant contains benzyl alcohol
                benzyl_alcohol_pattern = Chem.MolFromSmarts("c[CH2][OH]")
                has_benzyl_alcohol = any(
                    r.HasSubstructMatch(benzyl_alcohol_pattern) for r in reactants
                )

                # Check if product contains benzyl ether
                benzyl_ether_pattern = Chem.MolFromSmarts("c[CH2][O]c")
                has_benzyl_ether = (
                    product.HasSubstructMatch(benzyl_ether_pattern) if product else False
                )

                # If this is a benzyl ether formation reaction
                if has_phenol and has_benzyl_alcohol and has_benzyl_ether:
                    print("Found benzyl ether formation")
                    benzyl_ether_formation_found = True

        # Recursively process children
        for child in node.get("children", []):
            dfs_traverse(child)

    # Start traversal from the root
    dfs_traverse(route)

    return benzyl_ether_formation_found
