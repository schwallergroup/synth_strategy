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
    This function detects a synthetic strategy involving O-alkylation of a phenol
    with an alkyl halide to form an aryl ether.
    """
    found_o_alkylation = False

    def dfs_traverse(node):
        nonlocal found_o_alkylation

        if node["type"] == "reaction":
            if "metadata" in node and "rsmi" in node["metadata"]:
                rsmi = node["metadata"]["rsmi"]
                reactants = rsmi.split(">")[0].split(".")
                product = rsmi.split(">")[-1]

                # Check for phenol and alkyl halide
                has_phenol = False
                has_alkyl_halide = False

                for reactant in reactants:
                    try:
                        mol = Chem.MolFromSmiles(reactant)
                        if mol:
                            phenol_pattern = Chem.MolFromSmarts("[c][OH]")
                            alkyl_halide_pattern = Chem.MolFromSmarts("[C][Br,Cl,I]")

                            if mol.HasSubstructMatch(phenol_pattern):
                                has_phenol = True
                            if mol.HasSubstructMatch(alkyl_halide_pattern):
                                has_alkyl_halide = True
                    except:
                        pass

                # Check if product has aryl ether
                try:
                    prod_mol = Chem.MolFromSmiles(product)
                    if prod_mol:
                        aryl_ether_pattern = Chem.MolFromSmarts("[c][O][C]")
                        if prod_mol.HasSubstructMatch(aryl_ether_pattern):
                            if has_phenol and has_alkyl_halide:
                                found_o_alkylation = True
                                print("Found phenol O-alkylation")
                except:
                    pass

        # Traverse children
        for child in node.get("children", []):
            dfs_traverse(child)

    # Start traversal
    dfs_traverse(route)

    return found_o_alkylation
