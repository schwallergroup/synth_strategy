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
    This function detects biaryl coupling reactions involving halogenated aromatics.
    """
    biaryl_coupling_found = False

    def dfs_traverse(node):
        nonlocal biaryl_coupling_found

        if node["type"] == "reaction" and "metadata" in node and "rsmi" in node["metadata"]:
            rsmi = node["metadata"]["rsmi"]
            reactants = rsmi.split(">")[0].split(".")
            product = rsmi.split(">")[-1]

            # Check if reactants contain halogenated aromatics
            has_aryl_halide = False
            for reactant in reactants:
                if reactant:
                    mol = Chem.MolFromSmiles(reactant)
                    if mol:
                        aryl_br = mol.HasSubstructMatch(Chem.MolFromSmarts("[c][Br]"))
                        aryl_i = mol.HasSubstructMatch(Chem.MolFromSmarts("[c][I]"))
                        if aryl_br or aryl_i:
                            has_aryl_halide = True

            # Check if product has a biaryl bond that wasn't in reactants
            if has_aryl_halide and product:
                prod_mol = Chem.MolFromSmiles(product)
                if prod_mol:
                    # Look for biaryl bonds in product
                    if prod_mol.HasSubstructMatch(Chem.MolFromSmarts("[c]!@[c]")):
                        # Check if this bond wasn't in reactants
                        biaryl_in_reactants = False
                        for reactant in reactants:
                            if reactant:
                                r_mol = Chem.MolFromSmiles(reactant)
                                if r_mol and r_mol.HasSubstructMatch(
                                    Chem.MolFromSmarts("[c]!@[c]")
                                ):
                                    biaryl_in_reactants = True

                        if not biaryl_in_reactants:
                            print("Biaryl coupling via halides detected")
                            biaryl_coupling_found = True

        # Continue traversing
        for child in node.get("children", []):
            dfs_traverse(child)

    # Start traversal
    dfs_traverse(route)
    return biaryl_coupling_found
