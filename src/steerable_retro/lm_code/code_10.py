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
    This function detects a synthetic strategy involving coupling between an alkyne
    and an aromatic amine.
    """
    found_coupling = False

    def dfs_traverse(node):
        nonlocal found_coupling

        if node["type"] == "reaction":
            if "metadata" in node and "rsmi" in node["metadata"]:
                rsmi = node["metadata"]["rsmi"]
                reactants = rsmi.split(">")[0].split(".")

                # Check for alkyne and aromatic amine in reactants
                has_alkyne = False
                has_aromatic_amine = False

                for reactant in reactants:
                    try:
                        mol = Chem.MolFromSmiles(reactant)
                        if mol:
                            alkyne_pattern = Chem.MolFromSmarts("[C]#[C]")
                            aromatic_amine_pattern = Chem.MolFromSmarts("[c][NH2]")

                            if mol.HasSubstructMatch(alkyne_pattern):
                                has_alkyne = True
                            if mol.HasSubstructMatch(aromatic_amine_pattern):
                                has_aromatic_amine = True
                    except:
                        pass

                if has_alkyne and has_aromatic_amine:
                    found_coupling = True
                    print("Found alkyne-amine coupling")

        # Traverse children
        for child in node.get("children", []):
            dfs_traverse(child)

    # Start traversal
    dfs_traverse(route)

    return found_coupling
