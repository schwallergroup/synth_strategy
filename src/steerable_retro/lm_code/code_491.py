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
    Detects a synthetic strategy involving early bromination of an aromatic ring.
    """
    found_early_bromination = False

    # SMARTS pattern for aryl bromide
    aryl_bromide = Chem.MolFromSmarts("c[Br]")

    def dfs_traverse(node, depth=0):
        nonlocal found_early_bromination

        if node["type"] == "reaction" and depth >= 2:  # Early in synthesis
            # Extract reactants and product
            rsmi = node["metadata"]["rsmi"]
            reactants_smiles = rsmi.split(">")[0].split(".")
            product_smiles = rsmi.split(">")[-1]

            # Convert to RDKit molecules
            reactant_mols = [Chem.MolFromSmiles(r) for r in reactants_smiles if r]
            product_mol = Chem.MolFromSmiles(product_smiles)

            if product_mol and all(r for r in reactant_mols):
                # Check if this is a bromination reaction
                product_has_bromide = product_mol.HasSubstructMatch(aryl_bromide)
                reactants_have_bromide = any(
                    r.HasSubstructMatch(aryl_bromide) for r in reactant_mols
                )

                if product_has_bromide and not reactants_have_bromide:
                    print(f"Found early bromination at depth {depth}")
                    found_early_bromination = True

        # Continue traversing
        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    # Start traversal
    dfs_traverse(route)

    return found_early_bromination
