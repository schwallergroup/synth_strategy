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
    This function detects if the synthetic route involves late-stage addition
    of a complex aromatic moiety.
    """
    complex_aromatic_addition = False

    def dfs_traverse(node, depth=0):
        nonlocal complex_aromatic_addition

        if node["type"] == "reaction" and depth <= 1:  # Focus on late-stage reactions
            # Extract reactants and products
            rsmi = node["metadata"]["rsmi"]
            reactants = rsmi.split(">")[0].split(".")
            product = rsmi.split(">")[-1]

            # Check for complex aromatic systems in reactants
            complex_aromatic_pattern = Chem.MolFromSmarts("c1ccc(c2ccccc2)cc1")  # Biphenyl-like
            halogenated_aromatic_pattern = Chem.MolFromSmarts(
                "c-[Cl,Br,F,I]"
            )  # Halogenated aromatic

            for reactant in reactants:
                reactant_mol = Chem.MolFromSmiles(reactant)
                if reactant_mol:
                    # Check if it has complex aromatic system or halogenated aromatic
                    if reactant_mol.HasSubstructMatch(
                        complex_aromatic_pattern
                    ) or reactant_mol.HasSubstructMatch(halogenated_aromatic_pattern):
                        # Count aromatic rings
                        ring_info = reactant_mol.GetRingInfo()
                        aromatic_ring_count = 0

                        for ring in ring_info.AtomRings():
                            ring_atoms = [reactant_mol.GetAtomWithIdx(idx) for idx in ring]
                            if all(atom.GetIsAromatic() for atom in ring_atoms):
                                aromatic_ring_count += 1

                        if aromatic_ring_count >= 2:  # Complex if it has 2+ aromatic rings
                            print(f"Detected complex aromatic addition at depth {depth}")
                            complex_aromatic_addition = True

        # Traverse children
        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    # Start traversal
    dfs_traverse(route)

    if complex_aromatic_addition:
        print("Late-stage complex aromatic addition strategy detected")
    return complex_aromatic_addition
