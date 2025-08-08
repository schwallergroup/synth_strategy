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
    Detects a synthetic strategy involving late-stage reduction of C=C bond.
    """
    late_reduction_found = False

    def dfs_traverse(node, depth=0):
        nonlocal late_reduction_found

        if node["type"] == "reaction" and depth <= 1:  # Late stage = low depth
            rsmi = node["metadata"]["rsmi"]
            reactants = rsmi.split(">")[0]
            product = rsmi.split(">")[-1]

            # Check for C=C in reactants
            reactant_mol = Chem.MolFromSmiles(reactants)
            if reactant_mol and reactant_mol.HasSubstructMatch(Chem.MolFromSmarts("[#6]=[#6]")):
                # Check if the C=C is reduced in the product
                product_mol = Chem.MolFromSmiles(product)
                if product_mol:
                    # Count C=C in reactants and products
                    reactant_cc_count = len(
                        reactant_mol.GetSubstructMatches(Chem.MolFromSmarts("[#6]=[#6]"))
                    )
                    product_cc_count = len(
                        product_mol.GetSubstructMatches(Chem.MolFromSmarts("[#6]=[#6]"))
                    )

                    if reactant_cc_count > product_cc_count:
                        late_reduction_found = True
                        print(f"Found late-stage reduction at depth {depth}")

        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    dfs_traverse(route)

    print(f"Late-stage reduction strategy detected: {late_reduction_found}")
    return late_reduction_found
