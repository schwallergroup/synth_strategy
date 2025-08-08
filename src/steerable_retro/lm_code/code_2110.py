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
    This function detects a strategy involving nitration of an aromatic ring
    followed by reduction to an amine.
    """
    # Track nitration and reduction steps
    has_nitration = False
    has_nitro_reduction = False
    nitration_depth = -1
    reduction_depth = -1

    def dfs_traverse(node, depth=0):
        nonlocal has_nitration, has_nitro_reduction, nitration_depth, reduction_depth

        if node["type"] == "reaction":
            # Extract reactants and product
            rsmi = node["metadata"]["rsmi"]
            reactants_smiles = rsmi.split(">")[0].split(".")
            product_smiles = rsmi.split(">")[-1]

            # Convert to RDKit molecules
            reactants = [Chem.MolFromSmiles(r) for r in reactants_smiles]
            product = Chem.MolFromSmiles(product_smiles)

            if product and all(r for r in reactants):
                # Check for nitration (adding NO2 to aromatic ring)
                nitro_pattern = Chem.MolFromSmarts("c[N+](=[O])[O-]")
                aromatic_pattern = Chem.MolFromSmarts("c")

                if product.HasSubstructMatch(nitro_pattern) and not any(
                    r.HasSubstructMatch(nitro_pattern) for r in reactants
                ):
                    if any(r.HasSubstructMatch(aromatic_pattern) for r in reactants):
                        print(f"Detected nitration at depth {depth}")
                        has_nitration = True
                        nitration_depth = depth

                # Check for nitro reduction
                nitro_pattern = Chem.MolFromSmarts("[N+](=[O])[O-]")
                amine_pattern = Chem.MolFromSmarts("[NH2]")

                if any(
                    r.HasSubstructMatch(nitro_pattern) for r in reactants
                ) and product.HasSubstructMatch(amine_pattern):
                    print(f"Detected nitro reduction at depth {depth}")
                    has_nitro_reduction = True
                    reduction_depth = depth

        # Traverse children
        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    # Start traversal
    dfs_traverse(route)

    # Check if nitration occurs before reduction in the synthetic direction
    # (which means higher depth in retrosynthetic direction)
    return has_nitration and has_nitro_reduction and nitration_depth > reduction_depth
