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
    This function detects the construction of a polyether chain through
    sequential fragment assembly.
    """
    polyether_pattern = Chem.MolFromSmarts("[#6][#8][#6][#6][#8][#6]")

    found_polyether = False
    sequential_assembly = False
    ether_formation_steps = 0

    def dfs_traverse(node):
        nonlocal found_polyether, sequential_assembly, ether_formation_steps

        if node["type"] == "mol":
            mol = Chem.MolFromSmiles(node["smiles"])
            if mol and mol.HasSubstructMatch(polyether_pattern):
                found_polyether = True

        elif node["type"] == "reaction":
            if "rsmi" in node["metadata"]:
                rsmi = node["metadata"]["rsmi"]
                reactants = rsmi.split(">")[0].split(".")
                product = rsmi.split(">")[-1]

                product_mol = Chem.MolFromSmiles(product)

                # Count C-O-C formations
                if product_mol:
                    for reactant in reactants:
                        reactant_mol = Chem.MolFromSmiles(reactant)
                        if reactant_mol:
                            # Check if this reaction forms a new ether linkage
                            if product_mol.HasSubstructMatch(
                                Chem.MolFromSmarts("[#6][#8][#6]")
                            ) and not reactant_mol.HasSubstructMatch(
                                Chem.MolFromSmarts("[#6][#8][#6]")
                            ):
                                ether_formation_steps += 1
                                print("Found ether formation step")

        for child in node.get("children", []):
            dfs_traverse(child)

    dfs_traverse(route)

    # Consider it sequential assembly if multiple ether formation steps are found
    if ether_formation_steps >= 2:
        sequential_assembly = True

    return found_polyether and sequential_assembly
