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
    This function detects a synthetic strategy involving SNAr reaction for C-N bond formation.
    """
    snar_detected = False

    def dfs_traverse(node, depth=0):
        nonlocal snar_detected

        if node["type"] == "reaction":
            rsmi = node["metadata"]["rsmi"]
            reactants_smiles = rsmi.split(">")[0].split(".")
            product_smiles = rsmi.split(">")[-1]

            reactants = [Chem.MolFromSmiles(r) for r in reactants_smiles]
            product = Chem.MolFromSmiles(product_smiles)

            # Check for SNAr pattern - fluoronitrobenzene + amine
            fluoro_nitro_pattern = Chem.MolFromSmarts("Fc1c([N+](=O)[O-])cccc1")
            amine_pattern = Chem.MolFromSmarts("[NH2]c1ccccc1")
            cn_bond_pattern = Chem.MolFromSmarts("c1ccccc1Nc1ccccc1")

            reactants_with_fluoro_nitro = any(
                r and r.HasSubstructMatch(fluoro_nitro_pattern) for r in reactants
            )
            reactants_with_amine = any(r and r.HasSubstructMatch(amine_pattern) for r in reactants)

            if (
                reactants_with_fluoro_nitro
                and reactants_with_amine
                and product
                and product.HasSubstructMatch(cn_bond_pattern)
            ):
                print(f"Detected SNAr reaction at depth {depth}")
                snar_detected = True

        # Traverse children
        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    # Start traversal from the root
    dfs_traverse(route)

    print(f"SNAr for C-N bond formation detected: {snar_detected}")
    return snar_detected
