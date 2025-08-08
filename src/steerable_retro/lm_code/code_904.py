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
    This function detects if the route involves the formation of a diaryl ether bond.
    """
    has_diaryl_ether_formation = False

    def dfs_traverse(node):
        nonlocal has_diaryl_ether_formation

        if node["type"] == "reaction" and "rsmi" in node.get("metadata", {}):
            rsmi = node["metadata"]["rsmi"]
            reactants = rsmi.split(">")[0].split(".")
            product = rsmi.split(">")[-1]

            diaryl_ether_pattern = Chem.MolFromSmarts("[c]-[O]-[c]")

            product_mol = Chem.MolFromSmiles(product)

            if product_mol and product_mol.HasSubstructMatch(diaryl_ether_pattern):
                # Check if diaryl ether is not in any reactant
                diaryl_ether_in_reactants = False
                for reactant in reactants:
                    reactant_mol = Chem.MolFromSmiles(reactant)
                    if reactant_mol and reactant_mol.HasSubstructMatch(diaryl_ether_pattern):
                        diaryl_ether_in_reactants = True
                        break

                if not diaryl_ether_in_reactants:
                    has_diaryl_ether_formation = True
                    print("Diaryl ether formation detected")

        for child in node.get("children", []):
            dfs_traverse(child)

    dfs_traverse(route)
    print(f"Diaryl ether formation: {has_diaryl_ether_formation}")
    return has_diaryl_ether_formation
