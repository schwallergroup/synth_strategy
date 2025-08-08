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
    Detects if the synthesis route involves multiple C-N bond formations,
    including nucleophilic substitutions and reductive aminations.
    """
    cn_formation_count = 0

    def dfs_traverse(node):
        nonlocal cn_formation_count

        if node["type"] == "reaction" and "metadata" in node and "rsmi" in node["metadata"]:
            rsmi = node["metadata"]["rsmi"]
            reactants = rsmi.split(">")[0].split(".")
            product = rsmi.split(">")[-1]

            try:
                prod_mol = Chem.MolFromSmiles(product)

                # Check for C-N bond formation
                cn_pattern = Chem.MolFromSmarts("[#6]-[#7]")
                if prod_mol and prod_mol.HasSubstructMatch(cn_pattern):
                    # Count C-N bonds in product
                    prod_cn_count = len(prod_mol.GetSubstructMatches(cn_pattern))

                    # Count C-N bonds in reactants
                    react_cn_count = 0
                    for reactant in reactants:
                        react_mol = Chem.MolFromSmiles(reactant)
                        if react_mol:
                            react_cn_count += len(react_mol.GetSubstructMatches(cn_pattern))

                    # If product has more C-N bonds than reactants combined, a C-N bond was formed
                    if prod_cn_count > react_cn_count:
                        cn_formation_count += 1
                        print(
                            f"Found C-N bond formation at depth {node.get('metadata', {}).get('ID', '')}"
                        )
            except Exception as e:
                print(f"Error in SMARTS matching: {e}")

        for child in node.get("children", []):
            dfs_traverse(child)

    dfs_traverse(route)
    return cn_formation_count >= 2  # Return True if at least 2 C-N bonds were formed
