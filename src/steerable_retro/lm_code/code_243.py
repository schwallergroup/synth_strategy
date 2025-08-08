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
    This function detects if the synthesis route involves multiple sequential amide couplings.
    """
    amide_couplings = 0

    def dfs_traverse(node):
        nonlocal amide_couplings

        if node["type"] == "reaction" and "metadata" in node and "rsmi" in node["metadata"]:
            rsmi = node["metadata"]["rsmi"]
            reactants = rsmi.split(">")[0].split(".")
            product = rsmi.split(">")[-1]

            # Check for amide formation
            try:
                product_mol = Chem.MolFromSmiles(product)
                amide_pattern = Chem.MolFromSmarts("[#6]-[#6](=[#8])-[#7]")

                if product_mol and product_mol.HasSubstructMatch(amide_pattern):
                    # Count amide bonds in product
                    product_amide_count = len(product_mol.GetSubstructMatches(amide_pattern))

                    # Count amide bonds in reactants
                    reactant_amide_count = 0
                    for reactant in reactants:
                        reactant_mol = Chem.MolFromSmiles(reactant)
                        if reactant_mol:
                            reactant_amide_count += len(
                                reactant_mol.GetSubstructMatches(amide_pattern)
                            )

                    # If product has more amide bonds than reactants, an amide coupling occurred
                    if product_amide_count > reactant_amide_count:
                        print("Amide coupling detected")
                        amide_couplings += 1
            except:
                pass

        # Traverse children
        for child in node.get("children", []):
            dfs_traverse(child)

    # Call dfs_traverse on the root node
    dfs_traverse(route)

    # Return True if at least 2 amide couplings are found
    return amide_couplings >= 2
