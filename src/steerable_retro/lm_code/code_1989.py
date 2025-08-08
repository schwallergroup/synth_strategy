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
    Detects if the synthesis includes N-methylation of a sulfonamide.
    """
    methylation_found = False

    def dfs_traverse(node):
        nonlocal methylation_found

        if node["type"] == "reaction":
            if "rsmi" in node.get("metadata", {}):
                rsmi = node["metadata"]["rsmi"]
                reactants = rsmi.split(">")[0].split(".")
                product = rsmi.split(">")[-1]

                # Check for secondary sulfonamide in reactants
                secondary_sulfonamide_found = False
                for reactant in reactants:
                    reactant_mol = Chem.MolFromSmiles(reactant)
                    if reactant_mol:
                        secondary_sulfonamide_pattern = Chem.MolFromSmarts(
                            "[#6]-[#7H]-[#16](=[#8])(=[#8])-[#6]"
                        )
                        if reactant_mol.HasSubstructMatch(secondary_sulfonamide_pattern):
                            secondary_sulfonamide_found = True
                            break

                # Check for N-methylated sulfonamide in product
                product_mol = Chem.MolFromSmiles(product)
                if product_mol and secondary_sulfonamide_found:
                    tertiary_sulfonamide_pattern = Chem.MolFromSmarts(
                        "[#6]-[#7]([#6])-[#16](=[#8])(=[#8])-[#6]"
                    )
                    if product_mol.HasSubstructMatch(tertiary_sulfonamide_pattern):
                        print("Found N-methylation of sulfonamide")
                        methylation_found = True

        # Continue traversing
        for child in node.get("children", []):
            dfs_traverse(child)

    # Start traversal
    dfs_traverse(route)

    return methylation_found
