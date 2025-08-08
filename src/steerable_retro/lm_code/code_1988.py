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
    Detects if the synthesis follows a linear strategy (no convergent steps).

    A linear synthesis has no steps where multiple significant building blocks
    are combined. Significant building blocks are identified as child nodes in
    the synthesis tree.
    """
    is_linear = True

    def dfs_traverse(node):
        nonlocal is_linear

        if node["type"] == "reaction":
            if "rsmi" in node.get("metadata", {}):
                rsmi = node["metadata"]["rsmi"]
                reactants = rsmi.split(">")[0].split(".")

                # Get child node SMILES for comparison with reactants
                child_smiles = []
                for child in node.get("children", []):
                    if child["type"] == "mol" and "smiles" in child:
                        child_smiles.append(child["smiles"])

                # Count how many reactants correspond to child nodes
                child_reactants = 0
                for reactant in reactants:
                    # Check if this reactant corresponds to a child node
                    mol = Chem.MolFromSmiles(reactant)
                    if mol:
                        for child_smi in child_smiles:
                            child_mol = Chem.MolFromSmiles(child_smi)
                            if child_mol and Chem.MolToSmiles(mol) == Chem.MolToSmiles(child_mol):
                                child_reactants += 1
                                break

                if child_reactants > 1:
                    print(f"Found convergent step with {child_reactants} significant reactants")
                    is_linear = False

        # Continue traversing
        for child in node.get("children", []):
            dfs_traverse(child)

    # Start traversal
    dfs_traverse(route)

    if is_linear:
        print("Synthesis follows a linear strategy")

    return is_linear
