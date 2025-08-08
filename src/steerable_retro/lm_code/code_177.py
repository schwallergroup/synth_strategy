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
    Detects if the synthesis route uses a chlorination followed by displacement sequence.
    Looks for a carbonyl → chloro transformation followed by a chloro → OR transformation.
    """
    chlorination_reactions = []
    displacement_reactions = []

    def dfs_traverse(node, depth=0):
        if node["type"] == "reaction":
            # Extract reactants and product
            rsmi = node["metadata"]["rsmi"]
            product = rsmi.split(">")[-1]
            reactants = rsmi.split(">")[0].split(".")

            # Check for chlorination: C=O → C-Cl
            carbonyl_pattern = Chem.MolFromSmarts("[#6]=[#8]")
            chloro_pattern = Chem.MolFromSmarts("[#6]-[Cl]")

            # Check for displacement: C-Cl → C-OR
            ether_pattern = Chem.MolFromSmarts("[#6]-[#8]-[#6]")

            product_mol = Chem.MolFromSmiles(product)
            if not product_mol:
                return

            # Check for chlorination
            has_carbonyl_in_reactants = False
            for reactant in reactants:
                reactant_mol = Chem.MolFromSmiles(reactant)
                if reactant_mol and reactant_mol.HasSubstructMatch(carbonyl_pattern):
                    has_carbonyl_in_reactants = True
                    break

            if has_carbonyl_in_reactants and product_mol.HasSubstructMatch(chloro_pattern):
                chlorination_reactions.append((depth, product))
                print(f"Found chlorination reaction at depth {depth}")

            # Check for displacement
            has_chloro_in_reactants = False
            for reactant in reactants:
                reactant_mol = Chem.MolFromSmiles(reactant)
                if reactant_mol and reactant_mol.HasSubstructMatch(chloro_pattern):
                    has_chloro_in_reactants = True
                    break

            if has_chloro_in_reactants and product_mol.HasSubstructMatch(ether_pattern):
                displacement_reactions.append((depth, product))
                print(f"Found displacement reaction at depth {depth}")

        # Traverse children
        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    dfs_traverse(route)

    # Check if there's a chlorination followed by displacement
    if chlorination_reactions and displacement_reactions:
        for chlor_depth, _ in chlorination_reactions:
            for disp_depth, _ in displacement_reactions:
                if (
                    disp_depth < chlor_depth
                ):  # In retrosynthesis, lower depth means later in synthesis
                    print("Found chlorination-displacement sequence")
                    return True

    return False
