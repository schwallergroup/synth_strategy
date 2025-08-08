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
    This function detects vinyl group introduction via cross-coupling reaction.
    """
    vinyl_addition_detected = False

    def dfs_traverse(node):
        nonlocal vinyl_addition_detected

        if node["type"] == "reaction" and "metadata" in node and "rsmi" in node["metadata"]:
            rsmi = node["metadata"]["rsmi"]
            reactants = rsmi.split(">")[0].split(".")
            product = rsmi.split(">")[-1]

            try:
                product_mol = Chem.MolFromSmiles(product)

                # Check for aryl bromide in reactants
                aryl_bromide_pattern = Chem.MolFromSmarts("c[Br]")
                # Check for vinyl group in product
                vinyl_pattern = Chem.MolFromSmarts("C=C")
                # Check for boron reagent in reactants (Suzuki-type)
                boron_reagent = any("B" in r for r in reactants)

                aryl_bromide_in_reactants = any(
                    Chem.MolFromSmiles(r).HasSubstructMatch(aryl_bromide_pattern)
                    for r in reactants
                    if Chem.MolFromSmiles(r)
                )
                vinyl_in_product = product_mol and product_mol.HasSubstructMatch(vinyl_pattern)

                if aryl_bromide_in_reactants and vinyl_in_product and boron_reagent:
                    print("Detected vinyl addition via cross-coupling")
                    vinyl_addition_detected = True
            except:
                pass

        # Traverse children
        for child in node.get("children", []):
            dfs_traverse(child)

    # Start traversal from the root
    dfs_traverse(route)
    return vinyl_addition_detected
