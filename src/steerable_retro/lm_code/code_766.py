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
    Detects if the route has a functional group interconversion cascade
    (ester→acid→alcohol→chloride→N-alkylation).
    """
    # Track if we've seen each transformation
    transformations = {
        "ester_to_acid": False,
        "acid_to_alcohol": False,
        "alcohol_to_chloride": False,
        "chloride_to_n_alkylation": False,
    }

    def dfs_traverse(node, depth=0):
        if node["type"] == "reaction":
            # Extract reactants and product
            rsmi = node["metadata"]["rsmi"]
            reactants_smiles = rsmi.split(">")[0]
            product_smiles = rsmi.split(">")[-1]

            reactants = [Chem.MolFromSmiles(r) for r in reactants_smiles.split(".")]
            product = Chem.MolFromSmiles(product_smiles)

            # Define patterns
            ester_pattern = Chem.MolFromSmarts("[#6]-[#8]-[C]=O")
            acid_pattern = Chem.MolFromSmarts("[OH]-[C]=O")
            alcohol_pattern = Chem.MolFromSmarts("[#6]-[OH]")
            chloride_pattern = Chem.MolFromSmarts("[#6]-[Cl]")
            n_alkylation_pattern = Chem.MolFromSmarts("[#6]-[#7]")

            # Check for ester to acid
            if any(
                [r.HasSubstructMatch(ester_pattern) for r in reactants if r]
            ) and product.HasSubstructMatch(acid_pattern):
                print(f"Ester to acid transformation detected at depth {depth}")
                transformations["ester_to_acid"] = True

            # Check for acid to alcohol
            if any(
                [r.HasSubstructMatch(acid_pattern) for r in reactants if r]
            ) and product.HasSubstructMatch(alcohol_pattern):
                print(f"Acid to alcohol transformation detected at depth {depth}")
                transformations["acid_to_alcohol"] = True

            # Check for alcohol to chloride
            if any(
                [r.HasSubstructMatch(alcohol_pattern) for r in reactants if r]
            ) and product.HasSubstructMatch(chloride_pattern):
                print(f"Alcohol to chloride transformation detected at depth {depth}")
                transformations["alcohol_to_chloride"] = True

            # Check for chloride to N-alkylation
            if any(
                [r.HasSubstructMatch(chloride_pattern) for r in reactants if r]
            ) and product.HasSubstructMatch(n_alkylation_pattern):
                print(f"Chloride to N-alkylation transformation detected at depth {depth}")
                transformations["chloride_to_n_alkylation"] = True

        # Traverse children
        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    # Start traversal
    dfs_traverse(route)

    # Check if we have the complete cascade
    return all(transformations.values())
