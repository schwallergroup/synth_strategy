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
    This function detects SNAr ether formation between aryl halide and phenol.
    """
    ether_formation_via_snar = False

    def dfs_traverse(node):
        nonlocal ether_formation_via_snar

        if node["type"] == "reaction":
            if "metadata" in node and "rsmi" in node["metadata"]:
                rsmi = node["metadata"]["rsmi"]
                reactants = rsmi.split(">")[0].split(".")
                product = rsmi.split(">")[-1]

                # Check for aryl halide and phenol in reactants
                aryl_halide = False
                phenol = False

                for reactant in reactants:
                    if reactant:
                        mol = Chem.MolFromSmiles(reactant)
                        if mol:
                            aryl_halide_pattern = Chem.MolFromSmarts("c[Br,Cl,I,F]")
                            phenol_pattern = Chem.MolFromSmarts("c[OH]")

                            if mol.HasSubstructMatch(aryl_halide_pattern):
                                aryl_halide = True
                            if mol.HasSubstructMatch(phenol_pattern):
                                phenol = True

                # Check for diaryl ether in product
                if product and aryl_halide and phenol:
                    mol = Chem.MolFromSmiles(product)
                    if mol:
                        diaryl_ether_pattern = Chem.MolFromSmarts("c[#8]c")
                        if mol.HasSubstructMatch(diaryl_ether_pattern):
                            ether_formation_via_snar = True
                            print("SNAr ether formation detected")

        # Traverse children
        for child in node.get("children", []):
            dfs_traverse(child)

    # Start traversal
    dfs_traverse(route)

    return ether_formation_via_snar
