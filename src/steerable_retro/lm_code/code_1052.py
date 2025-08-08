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
    Detects if the synthetic route uses phenol O-alkylation to construct a linker
    between fragments (ArOH + X-R → ArOR)
    """
    found_o_alkylation = False

    def dfs_traverse(node, depth=0):
        nonlocal found_o_alkylation

        if node["type"] == "reaction":
            if "rsmi" in node["metadata"]:
                rsmi = node["metadata"]["rsmi"]
                reactants = rsmi.split(">")[0].split(".")
                product = rsmi.split(">")[-1]

                # Look for phenol pattern in reactants
                phenol_pattern = Chem.MolFromSmarts("c[OH]")
                alkyl_halide_pattern = Chem.MolFromSmarts("[CX4][Br,Cl,I]")

                has_phenol = False
                has_alkyl_halide = False

                for reactant in reactants:
                    mol = Chem.MolFromSmiles(reactant)
                    if mol:
                        if mol.HasSubstructMatch(phenol_pattern):
                            has_phenol = True
                        if mol.HasSubstructMatch(alkyl_halide_pattern):
                            has_alkyl_halide = True

                # Check if product has a new C-O bond
                if has_phenol and has_alkyl_halide:
                    # This is a simplified check - in a real implementation,
                    # you would need to track atom mappings to confirm the bond formation
                    found_o_alkylation = True
                    print(f"Found phenol O-alkylation at depth {depth}")

        # Traverse children
        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    # Start traversal
    dfs_traverse(route)

    return found_o_alkylation
