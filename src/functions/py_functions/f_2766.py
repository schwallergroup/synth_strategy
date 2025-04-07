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


def main(route):
    """
    This function detects synthesis of biaryl pyrazole compounds.
    """
    has_biaryl_pyrazole = False

    def dfs_traverse(node):
        nonlocal has_biaryl_pyrazole

        if node["type"] == "mol" and node.get("in_stock", False) == False:
            try:
                mol = Chem.MolFromSmiles(node["smiles"])
                if mol is None:
                    return

                # Check for pyrazole
                pyrazole_pattern = Chem.MolFromSmarts("[n]1[n]cc[c]1")

                # Check for biaryl system
                biaryl_pattern = Chem.MolFromSmarts("[c]-!@[c]")

                if mol.HasSubstructMatch(pyrazole_pattern) and mol.HasSubstructMatch(
                    biaryl_pattern
                ):
                    # Further check if one of the rings in the biaryl is the pyrazole
                    for bond in mol.GetBonds():
                        if bond.GetBondType() == rdkit.Chem.rdchem.BondType.SINGLE:
                            a1 = bond.GetBeginAtomIdx()
                            a2 = bond.GetEndAtomIdx()

                            # Check if both atoms are aromatic
                            if (
                                mol.GetAtomWithIdx(a1).GetIsAromatic()
                                and mol.GetAtomWithIdx(a2).GetIsAromatic()
                            ):
                                # Create fragments by breaking this bond
                                fragmentation = Chem.FragmentOnBonds(
                                    mol, [bond.GetIdx()], addDummies=False
                                )
                                frags = Chem.GetMolFrags(fragmentation, asMols=True)

                                if len(frags) == 2:
                                    # Check if one fragment has pyrazole
                                    if frags[0].HasSubstructMatch(
                                        pyrazole_pattern
                                    ) or frags[1].HasSubstructMatch(pyrazole_pattern):
                                        has_biaryl_pyrazole = True
                                        print("Detected biaryl pyrazole structure")
                                        break
            except Exception as e:
                print(f"Error processing molecule: {e}")

        # Process children
        for child in node.get("children", []):
            dfs_traverse(child)

    # Start traversal
    dfs_traverse(route)

    print(f"Biaryl pyrazole detected: {has_biaryl_pyrazole}")
    return has_biaryl_pyrazole
