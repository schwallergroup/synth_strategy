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
    This function detects a strategy where a morpholine-like ring is opened and later
    a new heterocyclic ring is formed.
    """
    # Track if we found the required patterns
    found_morpholine = False
    found_morpholine_opening = False
    found_new_heterocycle = False

    # SMARTS patterns
    morpholine_pattern = "[#7]1[#6][#6][#8][#6][#6]1"

    def dfs_traverse(node):
        nonlocal found_morpholine, found_morpholine_opening, found_new_heterocycle

        if node["type"] == "reaction":
            rsmi = node["metadata"]["rsmi"]
            reactants = rsmi.split(">")[0].split(".")
            product = rsmi.split(">")[-1]

            # Check for morpholine pattern
            if any(has_substructure(r, morpholine_pattern) for r in reactants):
                found_morpholine = True

                # Check if this is a ring opening reaction
                reactant_rings = sum(count_rings(r) for r in reactants)
                product_rings = count_rings(product)

                if reactant_rings > product_rings:
                    print(f"Found morpholine ring opening: {rsmi}")
                    found_morpholine_opening = True

            # Check for heterocycle formation
            if not any(has_substructure(r, morpholine_pattern) for r in reactants):
                reactant_rings = sum(count_rings(r) for r in reactants)
                product_rings = count_rings(product)

                if product_rings > reactant_rings:
                    # Check if product has a heterocycle
                    mol = Chem.MolFromSmiles(product)
                    if mol:
                        rings = Chem.GetSSSR(mol)
                        for ring in rings:
                            ring_atoms = [mol.GetAtomWithIdx(idx) for idx in ring]
                            if any(atom.GetSymbol() in ["N", "O", "S"] for atom in ring_atoms):
                                print(f"Found new heterocycle formation: {rsmi}")
                                found_new_heterocycle = True
                                break

        # Recursively process children
        for child in node.get("children", []):
            dfs_traverse(child)

    def has_substructure(smiles, pattern):
        mol = Chem.MolFromSmiles(smiles)
        if mol:
            patt = Chem.MolFromSmarts(pattern)
            return mol.HasSubstructMatch(patt)
        return False

    def count_rings(smiles):
        mol = Chem.MolFromSmiles(smiles)
        if mol:
            return len(Chem.GetSSSR(mol))
        return 0

    # Start traversal from the root
    dfs_traverse(route)

    # Return True if we found all required transformations
    return found_morpholine and found_morpholine_opening and found_new_heterocycle
