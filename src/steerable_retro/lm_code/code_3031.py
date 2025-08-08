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
    This function detects coupling between a heterocycle and piperidine fragment.
    """
    found_coupling = False

    def dfs_traverse(node):
        nonlocal found_coupling

        if node["type"] == "reaction" and "metadata" in node and "rsmi" in node["metadata"]:
            rsmi = node["metadata"]["rsmi"]
            reactants = rsmi.split(">")[0].split(".")
            product = rsmi.split(">")[-1]

            # Check for disconnection between heterocycle and piperidine
            if len(reactants) >= 2:  # Multiple fragments
                try:
                    product_mol = Chem.MolFromSmiles(product)

                    # Pattern for heterocycle-CH2-piperidine connection
                    coupling_pattern = Chem.MolFromSmarts("[c][CH2][N]1[CH2][CH2][CH][CH2][CH2]1")

                    # Check if product has the coupling pattern
                    if product_mol.HasSubstructMatch(coupling_pattern):
                        # Check if reactants have separate heterocycle and piperidine
                        piperidine_pattern = Chem.MolFromSmarts("[N]1[CH2][CH2][CH][CH2][CH2]1")
                        pyrazole_pattern = Chem.MolFromSmarts("[c]1[n][c][c][n][c]1")

                        found_piperidine = False
                        found_pyrazole = False

                        for reactant in reactants:
                            reactant_mol = Chem.MolFromSmiles(reactant)
                            if reactant_mol.HasSubstructMatch(piperidine_pattern):
                                found_piperidine = True
                            if reactant_mol.HasSubstructMatch(pyrazole_pattern):
                                found_pyrazole = True

                        if found_piperidine and found_pyrazole:
                            found_coupling = True
                            print("Found heterocycle-piperidine coupling")
                except:
                    pass

        for child in node.get("children", []):
            dfs_traverse(child)

    dfs_traverse(route)
    return found_coupling
