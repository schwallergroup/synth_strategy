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
    Detects if the synthesis route involves the formation of a pyrimidine ring.
    """
    pyrimidine_formation = False

    def dfs_traverse(node):
        nonlocal pyrimidine_formation

        if node["type"] == "reaction":
            if "rsmi" in node.get("metadata", {}):
                rsmi = node["metadata"]["rsmi"]
                reactants = rsmi.split(">")[0].split(".")
                product = rsmi.split(">")[-1]

                # Check if product contains pyrimidine but reactants don't
                product_mol = Chem.MolFromSmiles(product)
                reactant_mols = [Chem.MolFromSmiles(r) for r in reactants if r]

                pyrimidine_pattern = Chem.MolFromSmarts("[#6]1[#7][#6][#7][#6][#6]1")

                if product_mol and product_mol.HasSubstructMatch(pyrimidine_pattern):
                    has_pyrimidine_in_reactants = False
                    for mol in reactant_mols:
                        if mol and mol.HasSubstructMatch(pyrimidine_pattern):
                            has_pyrimidine_in_reactants = True
                            break

                    if not has_pyrimidine_in_reactants:
                        print("Detected pyrimidine ring formation")
                        pyrimidine_formation = True

        # Continue traversing
        for child in node.get("children", []):
            dfs_traverse(child)

    dfs_traverse(route)
    return pyrimidine_formation
