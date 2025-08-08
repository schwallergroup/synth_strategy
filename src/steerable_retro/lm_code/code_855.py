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
    This function detects if the synthetic route involves amide formation
    in the late stages of the synthesis.
    """
    amide_formation_depths = []

    def dfs_traverse(node, depth=0):
        if node["type"] == "reaction":
            if "rsmi" in node.get("metadata", {}):
                rsmi = node["metadata"]["rsmi"]
                reactants = rsmi.split(">")[0].split(".")
                product = rsmi.split(">")[-1]

                # Check for amide formation patterns
                product_mol = Chem.MolFromSmiles(product)
                amide_pattern = Chem.MolFromSmarts("C(=O)N")

                if product_mol and product_mol.HasSubstructMatch(amide_pattern):
                    # Check if this is a new amide formation
                    amine_pattern = Chem.MolFromSmarts("[NH2]")
                    acid_pattern = Chem.MolFromSmarts("C(=O)O")
                    acid_chloride_pattern = Chem.MolFromSmarts("C(=O)Cl")

                    has_amine = False
                    has_acid_or_derivative = False

                    for reactant in reactants:
                        reactant_mol = Chem.MolFromSmiles(reactant)
                        if reactant_mol:
                            if reactant_mol.HasSubstructMatch(amine_pattern):
                                has_amine = True
                            if reactant_mol.HasSubstructMatch(
                                acid_pattern
                            ) or reactant_mol.HasSubstructMatch(acid_chloride_pattern):
                                has_acid_or_derivative = True

                    if has_amine and has_acid_or_derivative:
                        amide_formation_depths.append(depth)
                        print(f"Amide formation detected at depth {depth}")

        # Traverse children
        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    # Start traversal
    dfs_traverse(route)

    # Check if there's an amide formation at depth 0-2 (late stage)
    return any(depth <= 2 for depth in amide_formation_depths)
