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
    This function detects if the route involves late-stage sulfonamide formation.
    Late-stage means in the final steps (low depth in retrosynthetic tree).
    """
    sulfonamide_formed = False
    late_stage = False

    def dfs_traverse(node, depth=0):
        nonlocal sulfonamide_formed, late_stage

        if node["type"] == "reaction" and depth <= 1:  # Late stage (final or penultimate step)
            if "rsmi" in node["metadata"]:
                rsmi = node["metadata"]["rsmi"]
                reactants = rsmi.split(">")[0].split(".")
                product = rsmi.split(">")[-1]

                # Check if product contains sulfonamide but reactants don't
                product_mol = Chem.MolFromSmiles(product)
                if product_mol and product_mol.HasSubstructMatch(
                    Chem.MolFromSmarts("[#7]-[#16](=[#8])(=[#8])-[#7]")
                ):
                    # Check if this is a formation (not present in reactants)
                    sulfonamide_in_reactants = False
                    for reactant in reactants:
                        reactant_mol = Chem.MolFromSmiles(reactant)
                        if reactant_mol and reactant_mol.HasSubstructMatch(
                            Chem.MolFromSmarts("[#7]-[#16](=[#8])(=[#8])-[#7]")
                        ):
                            sulfonamide_in_reactants = True
                            break

                    if not sulfonamide_in_reactants:
                        sulfonamide_formed = True
                        late_stage = True
                        print(f"Late-stage sulfonamide formation detected at depth {depth}")

        # Traverse children
        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    # Start traversal
    dfs_traverse(route)
    return sulfonamide_formed and late_stage
