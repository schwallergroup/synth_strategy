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
    This function detects if the synthetic route involves pyrazole ring formation.
    """
    pyrazole_formed = False

    def dfs_traverse(node):
        nonlocal pyrazole_formed

        if node["type"] == "reaction":
            if "metadata" in node and "rsmi" in node["metadata"]:
                rsmi = node["metadata"]["rsmi"]
                reactants_smiles = rsmi.split(">")[0]
                product_smiles = rsmi.split(">")[-1]

                # Check if pyrazole is in product but not in reactants
                pyrazole_pattern = "[#6]1[#7][#7][#6][#6]1"

                product_mol = Chem.MolFromSmiles(product_smiles)
                reactants_mol = Chem.MolFromSmiles(reactants_smiles)

                if product_mol and reactants_mol:
                    product_has_pyrazole = product_mol.HasSubstructMatch(
                        Chem.MolFromSmarts(pyrazole_pattern)
                    )
                    reactants_have_pyrazole = reactants_mol.HasSubstructMatch(
                        Chem.MolFromSmarts(pyrazole_pattern)
                    )

                    if product_has_pyrazole and not reactants_have_pyrazole:
                        print("Detected pyrazole formation")
                        pyrazole_formed = True

        for child in node.get("children", []):
            dfs_traverse(child)

    dfs_traverse(route)
    return pyrazole_formed
