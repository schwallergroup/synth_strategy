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
    This function detects a synthetic strategy involving multiple amide bond formations.
    """
    amide_coupling_count = 0

    def dfs_traverse(node):
        nonlocal amide_coupling_count

        if node["type"] == "reaction" and "metadata" in node and "rsmi" in node["metadata"]:
            rsmi = node["metadata"]["rsmi"]
            reactants = rsmi.split(">")[0].split(".")
            product = rsmi.split(">")[-1]

            # Check for carboxylic acid in reactants
            acid_pattern = Chem.MolFromSmarts("[#6](=[#8])[#8;H1]")

            # Check for amine in reactants
            amine_pattern = Chem.MolFromSmarts("[#7;H1,H2]")

            # Check for amide in product
            amide_pattern = Chem.MolFromSmarts("[#6](=[#8])[#7]")

            has_acid = False
            has_amine = False
            for reactant in reactants:
                mol = Chem.MolFromSmiles(reactant)
                if mol:
                    if mol.HasSubstructMatch(acid_pattern):
                        has_acid = True
                    if mol.HasSubstructMatch(amine_pattern):
                        has_amine = True

            product_mol = Chem.MolFromSmiles(product)
            has_amide = False
            if product_mol:
                if product_mol.HasSubstructMatch(amide_pattern):
                    has_amide = True

            # If reactants have acid and amine, and product has amide, it's an amide coupling
            if has_acid and has_amine and has_amide:
                amide_coupling_count += 1
                print(f"Detected amide coupling: {rsmi}")

        # Traverse children
        for child in node.get("children", []):
            dfs_traverse(child)

    # Start traversal
    dfs_traverse(route)

    # Strategy is present if there are multiple amide couplings
    result = amide_coupling_count >= 2
    print(f"Multiple amide coupling strategy detected: {result} (Count: {amide_coupling_count})")
    return result
