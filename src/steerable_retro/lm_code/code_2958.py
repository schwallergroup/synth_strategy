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
    This function detects if a synthesis route maintains stereochemistry
    in a beta-lactam core throughout the synthesis.
    """
    # Track if we've seen a beta-lactam with stereocenters
    found_stereo_beta_lactam = False
    stereo_preserved = True

    # Beta-lactam pattern
    beta_lactam_pattern = Chem.MolFromSmarts("[#6]1[#7][#6](=[#8])[#6]1")

    def dfs_traverse(node):
        nonlocal found_stereo_beta_lactam, stereo_preserved

        if node["type"] == "mol":
            mol = Chem.MolFromSmiles(node["smiles"])
            if mol is None:
                return

            # Check if molecule contains beta-lactam
            if mol.HasSubstructMatch(beta_lactam_pattern):
                # Count stereocenters in the beta-lactam
                matches = mol.GetSubstructMatches(beta_lactam_pattern)
                for match in matches:
                    # Check atoms in the beta-lactam for stereochemistry
                    stereo_atoms = 0
                    for atom_idx in match:
                        atom = mol.GetAtomWithIdx(atom_idx)
                        if atom.GetChiralTag() != Chem.ChiralType.CHI_UNSPECIFIED:
                            stereo_atoms += 1

                    if stereo_atoms > 0:
                        if not found_stereo_beta_lactam:
                            found_stereo_beta_lactam = True
                            print(f"Found beta-lactam with {stereo_atoms} stereocenters")
                        else:
                            # We've seen a stereo beta-lactam before, check if count matches
                            if stereo_atoms == 0:
                                stereo_preserved = False
                                print("Stereochemistry not preserved in beta-lactam")

        # Process children
        for child in node.get("children", []):
            dfs_traverse(child)

    # Start traversal
    dfs_traverse(route)

    # Return true if we found a stereo beta-lactam and stereochemistry was preserved
    if found_stereo_beta_lactam and stereo_preserved:
        print("Beta-lactam stereochemistry preservation strategy detected")
        return True

    return False
