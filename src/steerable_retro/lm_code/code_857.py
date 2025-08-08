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
    This function detects if the synthetic route uses acid chloride formation
    followed by amide coupling.
    """
    acid_chloride_depths = []
    amide_formation_depths = []

    def dfs_traverse(node, depth=0):
        if node["type"] == "reaction":
            if "rsmi" in node.get("metadata", {}):
                rsmi = node["metadata"]["rsmi"]
                reactants = rsmi.split(">")[0].split(".")
                product = rsmi.split(">")[-1]

                # Check for acid chloride formation
                product_mol = Chem.MolFromSmiles(product)
                acid_chloride_pattern = Chem.MolFromSmarts("C(=O)Cl")
                acid_pattern = Chem.MolFromSmarts("C(=O)O")

                if product_mol and product_mol.HasSubstructMatch(acid_chloride_pattern):
                    for reactant in reactants:
                        reactant_mol = Chem.MolFromSmiles(reactant)
                        if reactant_mol and reactant_mol.HasSubstructMatch(acid_pattern):
                            acid_chloride_depths.append(depth)
                            print(f"Acid chloride formation detected at depth {depth}")
                            break

                # Check for amide formation from acid chloride
                amide_pattern = Chem.MolFromSmarts("C(=O)N")
                amine_pattern = Chem.MolFromSmarts("[NH2]")

                if product_mol and product_mol.HasSubstructMatch(amide_pattern):
                    has_acid_chloride = False
                    has_amine = False

                    for reactant in reactants:
                        reactant_mol = Chem.MolFromSmiles(reactant)
                        if reactant_mol:
                            if reactant_mol.HasSubstructMatch(acid_chloride_pattern):
                                has_acid_chloride = True
                            if reactant_mol.HasSubstructMatch(amine_pattern):
                                has_amine = True

                    if has_acid_chloride and has_amine:
                        amide_formation_depths.append(depth)
                        print(f"Amide formation from acid chloride detected at depth {depth}")

        # Traverse children
        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    # Start traversal
    dfs_traverse(route)

    # Check if there's an acid chloride formation followed by amide formation
    # This means there should be an acid chloride formation at a higher depth
    # followed by an amide formation at a lower depth
    for acid_depth in acid_chloride_depths:
        for amide_depth in amide_formation_depths:
            if amide_depth < acid_depth:  # Remember: lower depth = later in synthesis
                return True

    return False
