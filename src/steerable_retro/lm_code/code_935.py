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
    This function detects a synthetic strategy involving the use of benzyl ethers
    as protecting groups for phenols.
    """
    benzyl_protection_count = 0

    def dfs_traverse(node):
        nonlocal benzyl_protection_count

        if node["type"] == "reaction" and "metadata" in node and "rsmi" in node["metadata"]:
            rsmi = node["metadata"]["rsmi"]
            reactants = rsmi.split(">")[0].split(".")
            product = rsmi.split(">")[-1]

            try:
                product_mol = Chem.MolFromSmiles(product)

                # Check for phenol in reactants
                has_phenol = False
                for reactant in reactants:
                    reactant_mol = Chem.MolFromSmiles(reactant)
                    if reactant_mol and reactant_mol.HasSubstructMatch(
                        Chem.MolFromSmarts("[c][OH]")
                    ):
                        has_phenol = True
                        break

                # Check for benzyl ether in product
                if (
                    product_mol
                    and product_mol.HasSubstructMatch(Chem.MolFromSmarts("[c][O][CH2][c]"))
                    and has_phenol
                ):
                    benzyl_protection_count += 1
                    print(f"Detected benzyl protection in reaction: {rsmi}")
            except:
                pass

        # Process children
        for child in node.get("children", []):
            dfs_traverse(child)

    # Start traversal
    dfs_traverse(route)

    # Return True if at least 1 benzyl protection is detected
    result = benzyl_protection_count >= 1
    print(f"Benzyl protecting group strategy detected: {result} (count: {benzyl_protection_count})")
    return result
