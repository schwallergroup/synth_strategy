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
    This function detects if the synthesis involves nitro reduction to amine.
    """
    nitro_reduction_detected = False

    def dfs_traverse(node):
        nonlocal nitro_reduction_detected

        if node["type"] == "reaction" and "metadata" in node and "rsmi" in node["metadata"]:
            rsmi = node["metadata"]["rsmi"]
            reactants = rsmi.split(">")[0].split(".")
            product = rsmi.split(">")[-1]

            # Check if reactant contains nitro group and product contains amine at same position
            nitro_pattern = Chem.MolFromSmarts("[#7+](=[#8])[#8-]")
            amine_pattern = Chem.MolFromSmarts("[NH2]")

            product_mol = Chem.MolFromSmiles(product) if product else None

            for reactant in reactants:
                if not reactant:
                    continue
                reactant_mol = Chem.MolFromSmiles(reactant)
                if (
                    reactant_mol
                    and reactant_mol.HasSubstructMatch(nitro_pattern)
                    and product_mol
                    and product_mol.HasSubstructMatch(amine_pattern)
                ):
                    # This is a simplified check - ideally we would verify the amine is at the same position as the nitro
                    print("Detected nitro reduction to amine")
                    nitro_reduction_detected = True

        for child in node.get("children", []):
            dfs_traverse(child)

    dfs_traverse(route)
    return nitro_reduction_detected
