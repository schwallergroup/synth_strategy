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
    This function detects the formation of an aromatic thioether (Ar-S-Ar) bond in the synthesis.
    """
    thioether_formation = False
    thioether_pattern = Chem.MolFromSmarts("c-[#16]-c")

    def dfs_traverse(node):
        nonlocal thioether_formation

        if node["type"] == "reaction":
            if "rsmi" in node.get("metadata", {}):
                rsmi = node["metadata"]["rsmi"]
                reactants = rsmi.split(">")[0]
                product = rsmi.split(">")[-1]

                reactants_mol = Chem.MolFromSmiles(reactants)
                product_mol = Chem.MolFromSmiles(product)

                # Check if thioether is being formed
                if (
                    product_mol
                    and product_mol.HasSubstructMatch(thioether_pattern)
                    and (
                        not reactants_mol or not reactants_mol.HasSubstructMatch(thioether_pattern)
                    )
                ):
                    thioether_formation = True
                    print(f"Aromatic thioether formation detected at reaction: {rsmi}")

        for child in node.get("children", []):
            dfs_traverse(child)

    dfs_traverse(route)
    print(f"Aromatic thioether formation present: {thioether_formation}")
    return thioether_formation
