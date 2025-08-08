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
    Detects if the synthesis includes a phenol O-alkylation step
    """
    found_o_alkylation = False

    def dfs_traverse(node):
        nonlocal found_o_alkylation

        if node.get("type") == "reaction" and "metadata" in node and "rsmi" in node["metadata"]:
            rsmi = node["metadata"]["rsmi"]
            reactants = rsmi.split(">")[0].split(".")
            product = rsmi.split(">")[-1]

            # Check for phenol O-alkylation
            phenol_pattern = Chem.MolFromSmarts("[c:1][O;H1:2]")
            aryl_ether_pattern = Chem.MolFromSmarts("[c:1][O;H0:2][C]")
            alkyl_halide_pattern = Chem.MolFromSmarts("[C][Br,Cl,I]")

            product_mol = Chem.MolFromSmiles(product)
            if product_mol and product_mol.HasSubstructMatch(aryl_ether_pattern):
                has_phenol = False
                has_alkyl_halide = False

                for reactant in reactants:
                    reactant_mol = Chem.MolFromSmiles(reactant)
                    if reactant_mol:
                        if reactant_mol.HasSubstructMatch(phenol_pattern):
                            has_phenol = True
                        if reactant_mol.HasSubstructMatch(alkyl_halide_pattern):
                            has_alkyl_halide = True

                if has_phenol and has_alkyl_halide:
                    found_o_alkylation = True
                    print("Detected phenol O-alkylation")

        # Traverse children
        for child in node.get("children", []):
            dfs_traverse(child)

    # Start traversal
    dfs_traverse(route)

    return found_o_alkylation
