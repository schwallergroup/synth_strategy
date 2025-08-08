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
    This function detects if the synthetic route uses late-stage amide bond formation
    as a key strategy (in the last 1-2 steps).
    """
    amide_formation_depths = []

    def dfs_traverse(node, depth=0):
        if node["type"] == "reaction":
            # Extract reactants and products
            rsmi = node["metadata"]["rsmi"]
            reactants = rsmi.split(">")[0].split(".")
            product = rsmi.split(">")[-1]

            # Check for amide formation
            # In retrosynthesis, amide bond breaking: product has amine and carboxylic acid
            # that were connected in the reactant

            # Patterns for amine, carboxylic acid, and amide
            amine_pattern = Chem.MolFromSmarts("[#7;H1,H2]")
            carboxylic_pattern = Chem.MolFromSmarts("[#6](=[O])[O]")
            amide_pattern = Chem.MolFromSmarts("[#6](=[O])[#7]")

            # Check if reactant has amide but products have amine and carboxylic acid
            reactant_mol = Chem.MolFromSmiles(
                product
            )  # Note: in retrosynthesis, product is the reactant

            if reactant_mol and reactant_mol.HasSubstructMatch(amide_pattern):
                # Check if any product has amine or carboxylic acid
                has_amine = False
                has_carboxylic = False

                for reactant in reactants:  # In retrosynthesis, reactants are the products
                    mol = Chem.MolFromSmiles(reactant)
                    if mol:
                        if mol.HasSubstructMatch(amine_pattern):
                            has_amine = True
                        if mol.HasSubstructMatch(carboxylic_pattern):
                            has_carboxylic = True

                if has_amine and has_carboxylic:
                    print(f"Detected amide bond formation at depth {depth}")
                    amide_formation_depths.append(depth)

        # Traverse children
        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    # Start traversal
    dfs_traverse(route)

    # Check if amide formation occurs in the late stage (depth 0-1)
    late_stage_amide = any(depth <= 1 for depth in amide_formation_depths)
    if late_stage_amide:
        print("Late-stage amide formation strategy detected")
    return late_stage_amide
