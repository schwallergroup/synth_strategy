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
    Detects if the synthesis route contains multiple nucleophilic aromatic substitutions
    (C-F to C-N or C-Cl to C-N transformations).
    """
    nas_count = 0

    def dfs_traverse(node, depth=0):
        nonlocal nas_count

        if node["type"] == "reaction":
            if "rsmi" in node.get("metadata", {}):
                rsmi = node["metadata"]["rsmi"]
                reactants_part = rsmi.split(">")[0]
                product_part = rsmi.split(">")[-1]

                # Check for patterns indicating nucleophilic aromatic substitution
                try:
                    reactants = reactants_part.split(".")

                    # Look for aryl halide in reactants
                    aryl_f_pattern = Chem.MolFromSmarts("[c][F]")
                    aryl_cl_pattern = Chem.MolFromSmarts("[c][Cl]")

                    # Look for nitrogen nucleophile
                    nitrogen_pattern = Chem.MolFromSmarts("[N]")

                    has_aryl_halide = False
                    has_nitrogen = False

                    for reactant in reactants:
                        mol = Chem.MolFromSmiles(reactant)
                        if mol:
                            if mol.HasSubstructMatch(aryl_f_pattern) or mol.HasSubstructMatch(
                                aryl_cl_pattern
                            ):
                                has_aryl_halide = True
                            if mol.HasSubstructMatch(nitrogen_pattern):
                                has_nitrogen = True

                    # Check if product has C-N bond where C-F or C-Cl was
                    if has_aryl_halide and has_nitrogen:
                        product_mol = Chem.MolFromSmiles(product_part)
                        if product_mol and product_mol.HasSubstructMatch(
                            Chem.MolFromSmarts("[c][N]")
                        ):
                            print(f"Found nucleophilic aromatic substitution at depth {depth}")
                            nas_count += 1
                except:
                    pass

        # Traverse children
        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    # Start traversal from the root
    dfs_traverse(route)
    return nas_count >= 2  # Return True if at least 2 NAS reactions are found
