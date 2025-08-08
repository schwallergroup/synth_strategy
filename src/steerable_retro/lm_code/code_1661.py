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
    This function detects indole synthesis via nitro reduction and cyclization.
    """
    has_nitro_reduction = False
    has_indole_formation = False

    def dfs_traverse(node):
        nonlocal has_nitro_reduction, has_indole_formation

        if node["type"] == "reaction":
            if "metadata" in node and "rsmi" in node["metadata"]:
                rsmi = node["metadata"]["rsmi"]
                reactants = rsmi.split(">")[0].split(".")
                product = rsmi.split(">")[-1]

                # Check for nitro reduction
                nitro_pattern = Chem.MolFromSmarts("[N+](=O)[O-]")
                amine_pattern = Chem.MolFromSmarts("[NH2]")

                for reactant in reactants:
                    reactant_mol = Chem.MolFromSmiles(reactant)
                    if reactant_mol and reactant_mol.HasSubstructMatch(nitro_pattern):
                        product_mol = Chem.MolFromSmiles(product)
                        if product_mol and product_mol.HasSubstructMatch(amine_pattern):
                            has_nitro_reduction = True
                            print("Detected nitro reduction")

                # Check for indole formation
                # Look for reactions where an aniline and a carbonyl fragment combine to form an indole
                product_mol = Chem.MolFromSmiles(product)
                if product_mol:
                    indole_pattern = Chem.MolFromSmarts("c1ccc2[nH]ccc2c1")
                    if product_mol.HasSubstructMatch(indole_pattern):
                        aniline_found = False
                        carbonyl_found = False

                        for reactant in reactants:
                            reactant_mol = Chem.MolFromSmiles(reactant)
                            if reactant_mol:
                                aniline_pattern = Chem.MolFromSmarts("c1ccc(N)cc1")
                                carbonyl_pattern = Chem.MolFromSmarts("[CX3](=O)[#6]")

                                if reactant_mol.HasSubstructMatch(aniline_pattern):
                                    aniline_found = True
                                if reactant_mol.HasSubstructMatch(carbonyl_pattern):
                                    carbonyl_found = True

                        if aniline_found and carbonyl_found:
                            has_indole_formation = True
                            print("Detected indole formation")

        for child in node.get("children", []):
            dfs_traverse(child)

    dfs_traverse(route)
    return has_nitro_reduction and has_indole_formation
