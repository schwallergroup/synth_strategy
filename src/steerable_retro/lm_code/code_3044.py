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
    This function detects a strategy involving hydrogenation of an alkene (C=C) to an alkane (C-C).
    """
    found_hydrogenation = False

    def dfs_traverse(node, depth=0):
        nonlocal found_hydrogenation

        if node["type"] == "reaction":
            if "rsmi" in node.get("metadata", {}):
                rsmi = node["metadata"]["rsmi"]
                reactants = rsmi.split(">")[0].split(".")
                product = rsmi.split(">")[-1]

                # Check for alkene pattern in reactants
                reactant_mol = Chem.MolFromSmiles(reactants[0])
                if reactant_mol:
                    alkene_pattern = Chem.MolFromSmarts("[#6]=[#6]")
                    if reactant_mol.HasSubstructMatch(alkene_pattern):
                        # Check if the alkene is reduced to alkane in the product
                        product_mol = Chem.MolFromSmiles(product)
                        if product_mol:
                            # Look for the absence of the alkene in the product
                            if not product_mol.HasSubstructMatch(alkene_pattern):
                                # This is a simplification - ideally we would check that the specific
                                # alkene was converted to an alkane, not just disappeared
                                found_hydrogenation = True
                                print(f"Found alkene hydrogenation at depth {depth}")

        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    dfs_traverse(route)
    return found_hydrogenation
