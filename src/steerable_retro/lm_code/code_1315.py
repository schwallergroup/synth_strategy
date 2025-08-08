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
    This function detects if the synthetic route involves formation of a triazole ring.
    """
    triazole_formation_found = False

    def dfs_traverse(node):
        nonlocal triazole_formation_found

        if node["type"] == "reaction":
            if "rsmi" in node.get("metadata", {}):
                rsmi = node["metadata"]["rsmi"]
                reactants = rsmi.split(">")[0].split(".")
                product = rsmi.split(">")[-1]

                # Check if product contains a triazole ring
                triazole_pattern = Chem.MolFromSmarts("c1nn[c]n1")  # Simplified triazole pattern

                try:
                    prod_mol = Chem.MolFromSmiles(product)
                    if prod_mol and prod_mol.HasSubstructMatch(triazole_pattern):
                        # Check if reactants don't have the triazole
                        has_triazole_in_reactants = False
                        for reactant in reactants:
                            try:
                                react_mol = Chem.MolFromSmiles(reactant)
                                if react_mol and react_mol.HasSubstructMatch(triazole_pattern):
                                    has_triazole_in_reactants = True
                                    break
                            except:
                                continue

                        if not has_triazole_in_reactants:
                            print("Triazole formation detected")
                            triazole_formation_found = True
                except:
                    pass

        for child in node.get("children", []):
            dfs_traverse(child)

    dfs_traverse(route)
    return triazole_formation_found
