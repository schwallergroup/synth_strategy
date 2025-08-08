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
    This function detects a late-stage urea formation strategy where a urea group
    is introduced in one of the final steps of the synthesis.
    """
    urea_formation_found = False

    def dfs_traverse(node, depth=0):
        nonlocal urea_formation_found

        if (
            node["type"] == "reaction" and depth <= 1
        ):  # Only consider reactions at depth 0 or 1 (late stage)
            if "metadata" in node and "rsmi" in node["metadata"]:
                rsmi = node["metadata"]["rsmi"]
                reactants = rsmi.split(">")[0].split(".")
                product = rsmi.split(">")[-1]

                # Check for urea formation (N-C(=O)-N pattern)
                urea_pattern = Chem.MolFromSmarts("[#7]-[#6](=[#8])-[#7]")

                prod_mol = Chem.MolFromSmiles(product)
                if prod_mol and prod_mol.HasSubstructMatch(urea_pattern):
                    # Check if reactants don't have urea
                    has_urea_in_reactants = False
                    for reactant in reactants:
                        react_mol = Chem.MolFromSmiles(reactant)
                        if react_mol and react_mol.HasSubstructMatch(urea_pattern):
                            has_urea_in_reactants = True
                            break

                    if not has_urea_in_reactants:
                        urea_formation_found = True
                        print(f"Late-stage urea formation detected at depth {depth}")

        # Traverse children with increased depth
        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    # Start traversal
    dfs_traverse(route)

    return urea_formation_found
