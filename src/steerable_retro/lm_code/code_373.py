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
    Detects synthesis strategy involving halogenated heterocycles as key intermediates.
    """
    halogenated_heterocycle_used = False

    def dfs_traverse(node, depth=0):
        nonlocal halogenated_heterocycle_used

        if node["type"] == "mol":
            if "smiles" in node:
                mol = Chem.MolFromSmiles(node["smiles"])
                if mol is not None:
                    # Check for halogenated heterocycle
                    if mol.HasSubstructMatch(Chem.MolFromSmarts("[n]")) and mol.HasSubstructMatch(
                        Chem.MolFromSmarts("[c][Cl,Br,I,F]")
                    ):
                        print(f"Found halogenated heterocycle at depth {depth}")
                        halogenated_heterocycle_used = True

        elif node["type"] == "reaction":
            if "rsmi" in node.get("metadata", {}):
                rsmi = node["metadata"]["rsmi"]
                reactants_smiles = rsmi.split(">")[0].split(".")

                for r_smiles in reactants_smiles:
                    mol = Chem.MolFromSmiles(r_smiles)
                    if mol is not None:
                        # Check for halogenated heterocycle
                        if mol.HasSubstructMatch(
                            Chem.MolFromSmarts("[n]")
                        ) and mol.HasSubstructMatch(Chem.MolFromSmarts("[c][Cl,Br,I,F]")):
                            print(f"Found halogenated heterocycle reactant at depth {depth}")
                            halogenated_heterocycle_used = True

        # Traverse children
        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    # Start traversal
    dfs_traverse(route)
    return halogenated_heterocycle_used
