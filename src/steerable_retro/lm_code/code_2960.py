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
    This function detects the use of an isocyanate intermediate for amide bond formation.
    """
    found_isocyanate_reaction = False

    def dfs_traverse(node):
        nonlocal found_isocyanate_reaction

        if node.get("type") == "reaction":
            if "metadata" in node and "rsmi" in node.get("metadata", {}):
                rsmi = node["metadata"]["rsmi"]
                reactants = rsmi.split(">")[0].split(".")

                # Check for isocyanate in reactants
                for r in reactants:
                    r_mol = Chem.MolFromSmiles(r)
                    if r_mol:
                        isocyanate_pattern = Chem.MolFromSmarts("[N]=[C]=[O]")
                        if r_mol.HasSubstructMatch(isocyanate_pattern):
                            print(
                                f"Detected isocyanate reactant at depth {node.get('depth', 'unknown')}"
                            )
                            found_isocyanate_reaction = True
                            break

        # Traverse children
        for child in node.get("children", []):
            dfs_traverse(child)

    # Start traversal
    dfs_traverse(route)

    print(f"Isocyanate-based amide formation: {found_isocyanate_reaction}")
    return found_isocyanate_reaction
