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
    Detects if the synthesis includes nucleophilic aromatic substitution steps,
    particularly chlorine displacement by an amine.
    """
    has_snar = False

    def dfs_traverse(node):
        nonlocal has_snar

        if node["type"] == "reaction":
            rsmi = node["metadata"].get("rsmi", "")
            reactants = rsmi.split(">")[0].split(".")
            products = rsmi.split(">")[-1]

            # Look for chloro-aryl in one reactant and amine in another
            has_chloroaryl = False
            has_amine = False

            for reactant in reactants:
                if "c" in reactant and "Cl" in reactant:
                    has_chloroaryl = True
                if "N" in reactant and not "Cl" in reactant:
                    has_amine = True

            # Check if product has C-N bond where chlorine was
            if (
                has_chloroaryl
                and has_amine
                and "c" in products
                and "N" in products
                and "Cl" not in products
            ):
                print("Found potential nucleophilic aromatic substitution")
                has_snar = True

        for child in node.get("children", []):
            dfs_traverse(child)

    dfs_traverse(route)
    return has_snar
