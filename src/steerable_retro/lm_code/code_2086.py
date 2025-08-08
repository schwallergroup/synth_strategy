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
    Detects if the synthesis includes a Williamson ether synthesis
    (formation of an aryl-alkyl ether bond from a phenol and an alkyl halide).
    """
    williamson_found = False

    def dfs_traverse(node, depth=0):
        nonlocal williamson_found

        if node["type"] == "reaction":
            if "metadata" in node and "rsmi" in node["metadata"]:
                rsmi = node["metadata"]["rsmi"]
                reactants = rsmi.split(">")[0].split(".")
                product = rsmi.split(">")[-1]

                # Check for Williamson ether synthesis
                product_mol = Chem.MolFromSmiles(product)
                aryl_alkyl_ether_pattern = Chem.MolFromSmarts("[#6;a]-[#8]-[#6;!a]")

                if product_mol and product_mol.HasSubstructMatch(aryl_alkyl_ether_pattern):
                    # Check if reactants include phenol and alkyl halide
                    phenol_pattern = Chem.MolFromSmarts("[#6;a]-[#8H]")
                    alkyl_halide_pattern = Chem.MolFromSmarts("[#6]-[Br,Cl,I,F]")

                    has_phenol = False
                    has_alkyl_halide = False

                    for reactant in reactants:
                        r_mol = Chem.MolFromSmiles(reactant)
                        if r_mol:
                            if r_mol.HasSubstructMatch(phenol_pattern):
                                has_phenol = True
                            if r_mol.HasSubstructMatch(alkyl_halide_pattern):
                                has_alkyl_halide = True

                    if has_phenol and has_alkyl_halide:
                        williamson_found = True
                        print(f"Williamson ether synthesis detected at depth {depth}")

        # Continue traversing the tree
        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    # Start traversal from the root
    dfs_traverse(route)

    return williamson_found
