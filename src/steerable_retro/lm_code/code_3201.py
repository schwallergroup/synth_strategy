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
    This function detects if the synthesis preserves stereochemistry throughout
    the route, maintaining chiral centers.
    """
    preserves_stereochemistry = True

    def dfs_traverse(node, depth=0):
        nonlocal preserves_stereochemistry

        if node["type"] == "reaction":
            if "rsmi" in node["metadata"]:
                rsmi = node["metadata"]["rsmi"]
                reactants = rsmi.split(">")[0].split(".")
                product = rsmi.split(">")[-1]

                # Count chiral centers in product
                product_mol = Chem.MolFromSmiles(product)
                if product_mol:
                    Chem.AssignStereochemistry(product_mol)
                    product_chiral_centers = len(Chem.FindMolChiralCenters(product_mol))

                    # Count chiral centers in reactants
                    reactant_chiral_centers = 0
                    for reactant in reactants:
                        reactant_mol = Chem.MolFromSmiles(reactant)
                        if reactant_mol:
                            Chem.AssignStereochemistry(reactant_mol)
                            reactant_chiral_centers += len(Chem.FindMolChiralCenters(reactant_mol))

                    # If product has fewer chiral centers than reactants combined,
                    # stereochemistry might not be preserved
                    if product_chiral_centers < reactant_chiral_centers:
                        print(f"Potential loss of stereochemistry at depth {depth}")
                        print(f"Product chiral centers: {product_chiral_centers}")
                        print(f"Reactant chiral centers: {reactant_chiral_centers}")
                        preserves_stereochemistry = False

        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    dfs_traverse(route)
    return preserves_stereochemistry
