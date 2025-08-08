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
    Detects reduction of nitro group to amine in the synthesis.
    """
    found_nitro_reduction = False

    def dfs_traverse(node):
        nonlocal found_nitro_reduction

        if node["type"] == "reaction" and node.get("metadata", {}).get("rsmi"):
            rsmi = node["metadata"]["rsmi"]
            reactants_smiles = rsmi.split(">")[0].split(".")
            product_smiles = rsmi.split(">")[-1]

            try:
                reactants = [Chem.MolFromSmiles(r) for r in reactants_smiles]
                product = Chem.MolFromSmiles(product_smiles)

                if product:
                    # Check for nitro pattern in reactants
                    nitro_pattern = Chem.MolFromSmarts("[#6][N+](=[O])[O-]")
                    has_nitro = any(
                        mol.HasSubstructMatch(nitro_pattern) for mol in reactants if mol
                    )

                    # Check for amine pattern in product
                    amine_pattern = Chem.MolFromSmarts("[#6][NH2]")
                    has_amine = product.HasSubstructMatch(amine_pattern)

                    if has_nitro and has_amine:
                        print("Found nitro to amine reduction")
                        found_nitro_reduction = True
            except:
                print("Error processing reaction SMILES for nitro reduction detection")

        # Process children
        for child in node.get("children", []):
            dfs_traverse(child)

    # Start traversal
    dfs_traverse(route)
    return found_nitro_reduction
