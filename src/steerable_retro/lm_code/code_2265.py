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
