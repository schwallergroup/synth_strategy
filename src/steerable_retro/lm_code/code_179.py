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
    Detects if the synthesis route involves multiple O-functionalization reactions
    (O-alkylation, O-methylation, etc.)
    """
    o_functionalization_count = 0

    def dfs_traverse(node):
        nonlocal o_functionalization_count

        if node["type"] == "reaction":
            # Extract reactants and product
            rsmi = node["metadata"]["rsmi"]
            product = rsmi.split(">")[-1]
            reactants = rsmi.split(">")[0].split(".")

            # Patterns for O-functionalization
            phenol_pattern = Chem.MolFromSmarts("[#8H]-[c]")
            ether_pattern = Chem.MolFromSmarts("[c]-[#8]-[C]")

            # Check for O-functionalization: ArOH → ArOR
            product_mol = Chem.MolFromSmiles(product)
            if not product_mol:
                return

            if product_mol.HasSubstructMatch(ether_pattern):
                for reactant in reactants:
                    reactant_mol = Chem.MolFromSmiles(reactant)
                    if reactant_mol and reactant_mol.HasSubstructMatch(phenol_pattern):
                        o_functionalization_count += 1
                        print(
                            f"Found O-functionalization reaction (count: {o_functionalization_count})"
                        )
                        break

        # Traverse children
        for child in node.get("children", []):
            dfs_traverse(child)

    dfs_traverse(route)
    return o_functionalization_count >= 2  # Return True if 2+ O-functionalization reactions
