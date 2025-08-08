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
    Detects if the synthesis involves multiple C-N bond formations.
    """
    cn_bond_formations = 0

    def dfs_traverse(node, depth=0):
        nonlocal cn_bond_formations

        if node["type"] == "reaction":
            if "rsmi" in node.get("metadata", {}):
                rsmi = node["metadata"]["rsmi"]
                reactants = rsmi.split(">")[0].split(".")
                product = rsmi.split(">")[-1]

                try:
                    # This is a simplified approach - in practice would need reaction mapping
                    product_mol = Chem.MolFromSmiles(product)
                    cn_bonds_product = 0
                    if product_mol:
                        cn_pattern = Chem.MolFromSmarts("[#6]-[#7]")
                        cn_bonds_product = len(product_mol.GetSubstructMatches(cn_pattern))

                    total_cn_bonds_reactants = 0
                    for reactant in reactants:
                        reactant_mol = Chem.MolFromSmiles(reactant)
                        if reactant_mol:
                            cn_pattern = Chem.MolFromSmarts("[#6]-[#7]")
                            total_cn_bonds_reactants += len(
                                reactant_mol.GetSubstructMatches(cn_pattern)
                            )

                    # If product has more C-N bonds than reactants combined, C-N bond formation occurred
                    if cn_bonds_product > total_cn_bonds_reactants:
                        cn_bond_formations += 1
                        print(f"Found C-N bond formation at depth {depth}")
                except:
                    pass

        # Traverse children
        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    dfs_traverse(route)
    return cn_bond_formations >= 3  # Return True if at least 3 C-N bonds are formed
