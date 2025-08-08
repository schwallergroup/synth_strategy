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
    This function detects a strategy involving multiple phenol O-functionalization steps
    (both O-methylation and O-alkylation).
    """
    # Track phenol functionalizations
    o_methylation = False
    o_alkylation = False

    def dfs_traverse(node):
        nonlocal o_methylation, o_alkylation

        if node["type"] == "reaction":
            # Extract reactants and product
            rsmi = node["metadata"]["rsmi"]
            reactants_smiles = rsmi.split(">")[0].split(".")
            product_smiles = rsmi.split(">")[-1]

            # Convert to RDKit molecules
            reactants = [Chem.MolFromSmiles(r) for r in reactants_smiles if r]
            product = Chem.MolFromSmiles(product_smiles)

            if product and all(r for r in reactants):
                # Check for phenol in reactants
                phenol_pattern = Chem.MolFromSmarts("[OH]c1ccccc1")

                # Check for O-methylation
                methoxy_pattern = Chem.MolFromSmarts("c1ccccc1OC")
                if product.HasSubstructMatch(methoxy_pattern):
                    if any(r.HasSubstructMatch(phenol_pattern) for r in reactants):
                        print("Detected O-methylation of phenol")
                        o_methylation = True

                # Check for O-alkylation (non-methyl)
                alkoxy_pattern = Chem.MolFromSmarts("c1ccccc1OCC")
                if product.HasSubstructMatch(alkoxy_pattern):
                    if any(r.HasSubstructMatch(phenol_pattern) for r in reactants):
                        print("Detected O-alkylation of phenol")
                        o_alkylation = True

        # Traverse children
        for child in node.get("children", []):
            dfs_traverse(child)

    # Start traversal
    dfs_traverse(route)

    # Check if both types of phenol functionalization are present
    strategy_present = o_methylation and o_alkylation
    print(f"Phenol functionalization strategy detection result: {strategy_present}")
    return strategy_present
