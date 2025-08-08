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
    Detects a synthesis featuring late-stage introduction of a heteroaromatic ring,
    typically via nucleophilic aromatic substitution in the final steps.
    """
    # Track late-stage heteroaryl introduction
    late_stage_heteroaryl_introduction = False

    # Pattern for heteroaromatic rings
    heteroaromatic_pattern = Chem.MolFromSmarts("[a;!c]1[a][a][a][a][a]1")

    def dfs_traverse(node, depth=0):
        nonlocal late_stage_heteroaryl_introduction

        if node["type"] == "reaction" and depth <= 1:  # Late stage (depth 0 or 1)
            if "rsmi" in node.get("metadata", {}):
                rsmi = node["metadata"]["rsmi"]
                reactants = rsmi.split(">")[0].split(".")
                product = rsmi.split(">")[-1]

                # Check if any reactant has a heteroaromatic ring
                reactant_has_heteroaromatic = False
                for reactant in reactants:
                    mol = Chem.MolFromSmiles(reactant)
                    if mol and mol.HasSubstructMatch(heteroaromatic_pattern):
                        reactant_has_heteroaromatic = True

                # Check if product has a heteroaromatic ring
                product_mol = Chem.MolFromSmiles(product)
                if product_mol and product_mol.HasSubstructMatch(heteroaromatic_pattern):
                    if reactant_has_heteroaromatic:
                        # Check if this is likely an SNAr reaction (halogen on heteroaromatic)
                        halogen_heteroaromatic = Chem.MolFromSmarts("[a;!c][a]([F,Cl,Br,I])")
                        for reactant in reactants:
                            mol = Chem.MolFromSmiles(reactant)
                            if mol and mol.HasSubstructMatch(halogen_heteroaromatic):
                                late_stage_heteroaryl_introduction = True
                                print(
                                    f"Late-stage heteroaryl introduction detected at depth {depth}"
                                )

        # Traverse children
        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    # Start traversal
    dfs_traverse(route)

    if late_stage_heteroaryl_introduction:
        print("Late-stage heteroaryl introduction strategy detected")
        return True
    return False
