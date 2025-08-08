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
    This function detects a synthesis route with multiple aromatic C-N bond formations.
    """
    aromatic_cn_formations = 0

    def dfs_traverse(node):
        nonlocal aromatic_cn_formations

        if node["type"] == "reaction" and "metadata" in node and "rsmi" in node["metadata"]:
            rsmi = node["metadata"]["rsmi"]
            reactants_str = rsmi.split(">")[0]
            product_str = rsmi.split(">")[-1]

            # Check for aromatic C-N bond formation
            try:
                product_mol = Chem.MolFromSmiles(product_str)
                if product_mol:
                    aromatic_cn_pattern = Chem.MolFromSmarts("[c][N]")
                    matches = product_mol.GetSubstructMatches(aromatic_cn_pattern)

                    # Check if these bonds exist in reactants
                    for reactant in reactants_str.split("."):
                        reactant_mol = Chem.MolFromSmiles(reactant)
                        if reactant_mol:
                            reactant_matches = reactant_mol.GetSubstructMatches(aromatic_cn_pattern)
                            # If there are more aromatic C-N bonds in product than in any reactant,
                            # it suggests a new aromatic C-N bond was formed
                            if len(matches) > len(reactant_matches):
                                aromatic_cn_formations += 1
                                print("Aromatic C-N bond formation detected")
                                break
            except:
                pass

        # Recursively traverse children
        for child in node.get("children", []):
            dfs_traverse(child)

    # Start traversal from the root
    dfs_traverse(route)

    # Check if the strategy involves multiple aromatic C-N bond formations
    return aromatic_cn_formations >= 2
