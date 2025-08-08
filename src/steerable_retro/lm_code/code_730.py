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
    This function detects the formation of a lactam ring (cyclic amide) in the synthesis.
    """
    has_lactam_formation = False

    def dfs_traverse(node):
        nonlocal has_lactam_formation

        if node["type"] == "reaction":
            rsmi = node["metadata"].get("rsmi", "")
            if not rsmi:
                return

            reactants = rsmi.split(">")[0].split(".")
            product = rsmi.split(">")[-1]

            # Check for lactam formation
            lactam_pattern = Chem.MolFromSmarts("[#6]1~[#6]~[#6]~[#6]~[#7]~[#6](=[#8])~1")
            acyclic_amide_pattern = Chem.MolFromSmarts("[#7]-[#6](=[#8])")

            product_mol = Chem.MolFromSmiles(product)
            if product_mol and product_mol.HasSubstructMatch(lactam_pattern):
                # Check if reactants don't have the lactam
                lactam_in_reactants = False
                for reactant in reactants:
                    mol = Chem.MolFromSmiles(reactant)
                    if mol and mol.HasSubstructMatch(lactam_pattern):
                        lactam_in_reactants = True
                        break

                if not lactam_in_reactants:
                    # Check if at least one reactant has an acyclic amide or amine
                    for reactant in reactants:
                        mol = Chem.MolFromSmiles(reactant)
                        if mol and (
                            mol.HasSubstructMatch(acyclic_amide_pattern)
                            or mol.HasSubstructMatch(Chem.MolFromSmarts("[#7]"))
                        ):
                            has_lactam_formation = True
                            print(f"Detected lactam formation: {rsmi}")
                            break

        # Traverse children
        for child in node.get("children", []):
            dfs_traverse(child)

    # Start traversal
    dfs_traverse(route)
    return has_lactam_formation
