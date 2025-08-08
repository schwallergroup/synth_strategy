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
    This function detects if the synthesis involves multiple transformations
    of nitrogen-containing functional groups.
    """
    n_transformations = 0

    def dfs_traverse(node, depth=0):
        nonlocal n_transformations

        if node["type"] == "reaction":
            if "rsmi" in node.get("metadata", {}):
                rsmi = node["metadata"]["rsmi"]
                reactants = rsmi.split(">")[0].split(".")
                product = rsmi.split(">")[-1]

                for reactant in reactants:
                    reactant_mol = Chem.MolFromSmiles(reactant)
                    product_mol = Chem.MolFromSmiles(product)

                    if reactant_mol and product_mol:
                        # Check for nitrogen-containing functional groups
                        n_patterns = [
                            Chem.MolFromSmarts("[NX3]"),  # Amine
                            Chem.MolFromSmarts("[NX3H2]"),  # Primary amine
                            Chem.MolFromSmarts("[NX3H]([#6])[#6]"),  # Secondary amine
                            Chem.MolFromSmarts("[NX3]([#6])([#6])[#6]"),  # Tertiary amine
                            Chem.MolFromSmarts("[N+](=O)[O-]"),  # Nitro
                        ]

                        reactant_n_groups = sum(
                            [1 for pattern in n_patterns if reactant_mol.HasSubstructMatch(pattern)]
                        )
                        product_n_groups = sum(
                            [1 for pattern in n_patterns if product_mol.HasSubstructMatch(pattern)]
                        )

                        if (
                            reactant_n_groups > 0
                            and product_n_groups > 0
                            and reactant_n_groups != product_n_groups
                        ):
                            print(
                                f"Nitrogen functional group transformation detected at depth {depth}"
                            )
                            n_transformations += 1

        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    dfs_traverse(route)
    return n_transformations >= 2  # At least 2 transformations
