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
    This function detects a synthetic strategy involving reduction of a nitro group to an amine.
    """
    nitro_reduction = False

    def dfs_traverse(node):
        nonlocal nitro_reduction

        if node["type"] == "reaction":
            if "metadata" in node and "rsmi" in node["metadata"]:
                rsmi = node["metadata"]["rsmi"]
                reactants = rsmi.split(">")[0].split(".")
                product = rsmi.split(">")[-1]

                # Check for nitro group in reactant
                nitro_pattern = Chem.MolFromSmarts("[N+](=O)[O-]")
                nitro_present = False
                try:
                    for reactant in reactants:
                        mol = Chem.MolFromSmiles(reactant)
                        if mol and mol.HasSubstructMatch(nitro_pattern):
                            nitro_present = True
                            print("Nitro group detected in reactant")
                            break
                except:
                    pass

                # Check for amine in product at same position
                if nitro_present:
                    amine_pattern = Chem.MolFromSmarts("N")
                    try:
                        prod_mol = Chem.MolFromSmiles(product)
                        if prod_mol and prod_mol.HasSubstructMatch(amine_pattern):
                            # This is a simplification - ideally we would check that the amine is at the same position
                            # as the nitro group, but that requires atom mapping
                            nitro_reduction = True
                            print("Potential nitro reduction to amine detected")
                    except:
                        pass

        # Traverse children
        for child in node.get("children", []):
            dfs_traverse(child)

    # Start traversal
    dfs_traverse(route)

    return nitro_reduction
