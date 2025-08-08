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
    Detects if the synthetic route involves a quaternary carbon center
    connected to a cyclohexyl group.
    """
    quaternary_carbon_pattern = Chem.MolFromSmarts("[#6]([#6])([#6])([#6])[#6]")
    cyclohexyl_pattern = Chem.MolFromSmarts("[#6]1[#6][#6][#6][#6][#6]1")

    found_quaternary_carbon_cyclohexyl = False

    def dfs_traverse(node):
        nonlocal found_quaternary_carbon_cyclohexyl

        if node["type"] == "mol" and "smiles" in node:
            try:
                mol = Chem.MolFromSmiles(node["smiles"])
                if mol:
                    # Check if molecule has both quaternary carbon and cyclohexyl
                    quat_matches = mol.GetSubstructMatches(quaternary_carbon_pattern)
                    cyclohexyl_matches = mol.GetSubstructMatches(cyclohexyl_pattern)

                    if quat_matches and cyclohexyl_matches:
                        # Check if any quaternary carbon is connected to cyclohexyl
                        for quat_match in quat_matches:
                            quat_atom = mol.GetAtomWithIdx(quat_match[0])
                            for neighbor in quat_atom.GetNeighbors():
                                for cyclohexyl_match in cyclohexyl_matches:
                                    if neighbor.GetIdx() in cyclohexyl_match:
                                        print("Found quaternary carbon connected to cyclohexyl")
                                        found_quaternary_carbon_cyclohexyl = True
            except:
                pass

        for child in node.get("children", []):
            dfs_traverse(child)

    dfs_traverse(route)
    return found_quaternary_carbon_cyclohexyl
