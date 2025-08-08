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
    This function detects if the synthetic route includes a heterocycle formation step,
    specifically focusing on imidazole ring formation.
    """
    result = False

    def dfs_traverse(node):
        nonlocal result

        if node["type"] == "reaction" and node.get("metadata", {}).get("rsmi"):
            rsmi = node["metadata"]["rsmi"]
            reactants_smiles = rsmi.split(">")[0].split(".")
            product_smiles = rsmi.split(">")[-1]

            # Count rings in reactants and product
            product_mol = Chem.MolFromSmiles(product_smiles)

            if product_mol:
                # Check specifically for imidazole formation
                imidazole_pattern = Chem.MolFromSmarts("[#6]1[#6][#7][#6][#7]1")
                if product_mol.HasSubstructMatch(imidazole_pattern):
                    # Check if imidazole is formed in this step
                    imidazole_in_reactants = False
                    for reactant in reactants_smiles:
                        reactant_mol = Chem.MolFromSmiles(reactant)
                        if reactant_mol and reactant_mol.HasSubstructMatch(imidazole_pattern):
                            imidazole_in_reactants = True
                            break

                    if not imidazole_in_reactants:
                        print("Detected imidazole ring formation")
                        result = True

        for child in node.get("children", []):
            dfs_traverse(child)

    dfs_traverse(route)
    return result
