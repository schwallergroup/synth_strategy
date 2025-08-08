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
    Detects if the synthesis route includes coupling of pyridine and dimethylpyrrole moieties
    """
    pyridine_pattern = Chem.MolFromSmarts("c1ncccc1")
    dimethylpyrrole_pattern = Chem.MolFromSmarts("Cc1ccc(C)n1")

    def dfs_traverse(node):
        if node["type"] == "reaction" and "metadata" in node and "rsmi" in node["metadata"]:
            rsmi = node["metadata"]["rsmi"]
            product = rsmi.split(">")[-1]

            try:
                p_mol = Chem.MolFromSmiles(product)
                if (
                    p_mol
                    and p_mol.HasSubstructMatch(pyridine_pattern)
                    and p_mol.HasSubstructMatch(dimethylpyrrole_pattern)
                ):

                    # Check if this is a coupling reaction by examining reactants
                    reactants = rsmi.split(">")[0].split(".")
                    pyridine_in_reactants = False
                    dimethylpyrrole_in_reactants = False

                    for reactant in reactants:
                        try:
                            r_mol = Chem.MolFromSmiles(reactant)
                            if r_mol:
                                if r_mol.HasSubstructMatch(pyridine_pattern):
                                    pyridine_in_reactants = True
                                if r_mol.HasSubstructMatch(dimethylpyrrole_pattern):
                                    dimethylpyrrole_in_reactants = True
                        except:
                            continue

                    if pyridine_in_reactants and dimethylpyrrole_in_reactants:
                        print("Pyridine-dimethylpyrrole coupling detected")
                        return True
            except:
                pass

        for child in node.get("children", []):
            if dfs_traverse(child):
                return True

        return False

    return dfs_traverse(route)
