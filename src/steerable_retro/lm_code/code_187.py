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
    This function detects if the synthesis includes a malonic ester alkylation step.
    """
    has_malonic_alkylation = False

    def dfs_traverse(node):
        nonlocal has_malonic_alkylation

        if node["type"] == "reaction":
            if "rsmi" in node.get("metadata", {}):
                rsmi = node["metadata"]["rsmi"]
                reactants = rsmi.split(">")[0].split(".")
                product = rsmi.split(">")[-1]

                # Pattern for malonic ester (diethyl malonate or similar)
                malonic_pattern = Chem.MolFromSmarts("C(=O)OC.C(=O)OC")

                # Pattern for alkyl halide
                alkyl_halide_pattern = Chem.MolFromSmarts("[#6][Cl,Br,I]")

                # Check reactants
                has_malonic = False
                has_alkyl_halide = False

                for reactant in reactants:
                    mol = Chem.MolFromSmiles(reactant)
                    if not mol:
                        continue

                    if mol.HasSubstructMatch(malonic_pattern):
                        has_malonic = True

                    if mol.HasSubstructMatch(alkyl_halide_pattern):
                        has_alkyl_halide = True

                if has_malonic and has_alkyl_halide:
                    has_malonic_alkylation = True
                    print("Found malonic ester alkylation")

        # Continue traversing
        for child in node.get("children", []):
            dfs_traverse(child)

    dfs_traverse(route)
    print(f"Malonic ester alkylation: {'present' if has_malonic_alkylation else 'absent'}")
    return has_malonic_alkylation
