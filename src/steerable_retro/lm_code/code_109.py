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
    This function detects Suzuki coupling reactions in the synthesis.
    """
    suzuki_coupling_found = False

    def dfs_traverse(node):
        nonlocal suzuki_coupling_found

        if node["type"] == "reaction":
            if "rsmi" in node.get("metadata", {}):
                rsmi = node["metadata"]["rsmi"]
                reactants = rsmi.split(">")[0].split(".")
                product = rsmi.split(">")[-1]

                # Check for boronic acid/ester in reactants
                boronic_acid_present = False
                aryl_halide_present = False

                for reactant in reactants:
                    mol = Chem.MolFromSmiles(reactant)
                    if mol:
                        # Check for boronic acid/ester
                        boronic_pattern = Chem.MolFromSmarts("[c]-[B]([OX2])[OX2]")
                        if mol.HasSubstructMatch(boronic_pattern):
                            boronic_acid_present = True

                        # Check for aryl halide/pseudohalide
                        aryl_halide_pattern = Chem.MolFromSmarts("[c]-[Br,I,Cl,F]")
                        if mol.HasSubstructMatch(aryl_halide_pattern):
                            aryl_halide_present = True

                # Check if product has a new biaryl bond
                if boronic_acid_present and aryl_halide_present:
                    prod_mol = Chem.MolFromSmiles(product)
                    if prod_mol:
                        biaryl_pattern = Chem.MolFromSmarts("c-c")
                        if prod_mol.HasSubstructMatch(biaryl_pattern):
                            suzuki_coupling_found = True

        # Traverse children
        for child in node.get("children", []):
            dfs_traverse(child)

    # Start traversal
    dfs_traverse(route)

    print(f"Suzuki coupling detected: {suzuki_coupling_found}")
    return suzuki_coupling_found
