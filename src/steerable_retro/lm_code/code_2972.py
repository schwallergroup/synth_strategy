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
    This function detects if an amide formation occurs in the early stages of the synthesis.
    """
    early_amide_formation = False

    def dfs_traverse(node, depth=0):
        nonlocal early_amide_formation

        if node["type"] == "reaction" and depth >= 3:  # Early stage (depth 3 or higher)
            if "rsmi" in node.get("metadata", {}):
                rsmi = node["metadata"]["rsmi"]
                reactants = rsmi.split(">")[0].split(".")
                product = rsmi.split(">")[-1]

                # Check for acyl chloride pattern in reactants
                acyl_chloride_pattern = Chem.MolFromSmarts("[#6](=[#8])[Cl]")
                amine_pattern = Chem.MolFromSmarts("[#7;H1,H2]")
                amide_pattern = Chem.MolFromSmarts("[#6](=[#8])[#7]")

                # Check if reactants contain acyl chloride and amine
                acyl_chloride_present = False
                amine_present = False

                for reactant in reactants:
                    reactant_mol = Chem.MolFromSmiles(reactant)
                    if reactant_mol:
                        if reactant_mol.HasSubstructMatch(acyl_chloride_pattern):
                            acyl_chloride_present = True
                        if reactant_mol.HasSubstructMatch(amine_pattern):
                            amine_present = True

                # Check if product contains amide
                product_mol = Chem.MolFromSmiles(product)
                amide_in_product = product_mol and product_mol.HasSubstructMatch(amide_pattern)

                if acyl_chloride_present and amine_present and amide_in_product:
                    early_amide_formation = True
                    print(f"Detected amide formation at depth {depth}")

        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    dfs_traverse(route)
    return early_amide_formation
