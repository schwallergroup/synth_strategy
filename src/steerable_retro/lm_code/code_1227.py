#!/bin/python

"""LM-defined function for strategy description."""

import copy
import re
from collections import deque

import rdkit
import rdkit.Chem as Chem
from rdkit import Chem
from rdkit.Chem import (
    AllChem,
    Descriptors,
    Lipinski,
    rdChemReactions,
    rdFMCS,
    rdMolDescriptors,
    rdmolops,
)
from rdkit.Chem.Scaffolds import MurckoScaffold


def main(smiles):
    """
    Determines if a molecule is likely an auxiliary reagent (solvent, catalyst, etc.)
    based on simple heuristics.
    """
    # Common simple molecules that are often auxiliary reagents
    simple_molecules = ["O", "I", "Cl", "Br", "F", "H2O", "HCl", "HBr", "HI", "NH3"]
    if smiles in simple_molecules:
        return True

    # Check for very small molecules (likely reagents, not main compounds)
    mol = Chem.MolFromSmiles(smiles)
    if mol and mol.GetNumAtoms() < 5:
        return True

    return False
