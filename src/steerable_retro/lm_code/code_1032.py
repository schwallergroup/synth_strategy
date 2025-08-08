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
    This function detects if the synthesis involves introduction of a halogen (Br, Cl, I, F)
    in the late stage (first 2 steps).
    """
    late_halogenation_found = False

    def dfs_traverse(node, depth=0):
        nonlocal late_halogenation_found

        if (
            node["type"] == "reaction" and depth <= 1
        ):  # Focus on very late-stage reactions (depth 0 or 1)
            if "metadata" in node and "rsmi" in node["metadata"]:
                rsmi = node["metadata"]["rsmi"]
                reactants_smiles = rsmi.split(">")[0]
                product_smiles = rsmi.split(">")[-1]

                reactants_mol = Chem.MolFromSmiles(reactants_smiles)
                product_mol = Chem.MolFromSmiles(product_smiles)

                if reactants_mol and product_mol:
                    # Check for halogen introduction
                    halogen_smarts = "[#6]-[#9,#17,#35,#53]"  # C-F, C-Cl, C-Br, C-I

                    reactants_has_halogen = reactants_mol.HasSubstructMatch(
                        Chem.MolFromSmarts(halogen_smarts)
                    )
                    product_has_halogen = product_mol.HasSubstructMatch(
                        Chem.MolFromSmarts(halogen_smarts)
                    )

                    if not reactants_has_halogen and product_has_halogen:
                        late_halogenation_found = True
                        print(f"Detected late-stage halogenation at depth {depth}")

        # Continue traversing
        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    # Start traversal from the root
    dfs_traverse(route)

    return late_halogenation_found
