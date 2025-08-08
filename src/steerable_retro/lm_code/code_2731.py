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
    This function detects if the synthesis involves a sequence of at least 2 functional group
    interconversions (e.g., ester → alcohol → chloride → nitrile).
    """
    # Define patterns for functional groups
    patterns = {
        "ester": Chem.MolFromSmarts("C(=O)OC"),
        "alcohol": Chem.MolFromSmarts("CO"),
        "chloride": Chem.MolFromSmarts("CCl"),
        "nitrile": Chem.MolFromSmarts("C#N"),
    }

    # Track functional group transformations
    transformations = []

    def dfs_traverse(node):
        if node["type"] == "reaction" and "rsmi" in node.get("metadata", {}):
            rsmi = node["metadata"]["rsmi"]
            reactants_part = rsmi.split(">")[0]
            product_part = rsmi.split(">")[-1]

            reactant_mol = Chem.MolFromSmiles(reactants_part)
            product_mol = Chem.MolFromSmiles(product_part)

            if reactant_mol and product_mol:
                # Check which functional groups are present in reactants and products
                reactant_groups = []
                product_groups = []

                for group_name, pattern in patterns.items():
                    if reactant_mol.HasSubstructMatch(pattern):
                        reactant_groups.append(group_name)
                    if product_mol.HasSubstructMatch(pattern):
                        product_groups.append(group_name)

                # If there's a change in functional groups, record it
                if set(reactant_groups) != set(product_groups):
                    transformations.append((reactant_groups, product_groups))

        for child in node.get("children", []):
            dfs_traverse(child)

    dfs_traverse(route)

    # Count distinct functional group transformations
    distinct_transformations = len(transformations)
    result = distinct_transformations >= 2

    print(
        f"Sequential functional group interconversions: {result} (detected {distinct_transformations} transformations)"
    )
    return result
