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
    Detects a strategy where heterocyclic cores are decorated with
    nitrogen-containing heterocycles via SNAr reactions.
    """
    # Initialize tracking variables
    heterocycle_decoration_count = 0
    heterocycles_present = set()

    def has_heterocycle(smiles):
        mol = Chem.MolFromSmiles(smiles)
        if not mol:
            return False

        # Check for common heterocycles
        thiazole_pattern = Chem.MolFromSmarts("c1sccn1")
        piperazine_pattern = Chem.MolFromSmarts("N1CCNCC1")
        morpholine_pattern = Chem.MolFromSmarts("N1CCOCC1")

        if mol.HasSubstructMatch(thiazole_pattern):
            heterocycles_present.add("thiazole")
            return True
        if mol.HasSubstructMatch(piperazine_pattern):
            heterocycles_present.add("piperazine")
            return True
        if mol.HasSubstructMatch(morpholine_pattern):
            heterocycles_present.add("morpholine")
            return True

        return False

    def is_heterocycle_decoration(reactants, product):
        # Check if a heterocycle is being attached to another ring system
        heterocycle_reactant = False
        aromatic_reactant = False

        for reactant in reactants:
            mol = Chem.MolFromSmiles(reactant)
            if not mol:
                continue

            if mol.HasSubstructMatch(Chem.MolFromSmarts("[a]")):  # Aromatic atom
                aromatic_reactant = True

            if has_heterocycle(reactant):
                heterocycle_reactant = True

        # Check if product has both heterocycle and aromatic system
        prod_mol = Chem.MolFromSmiles(product)
        if heterocycle_reactant and aromatic_reactant and prod_mol:
            # Check if product has a new C-N bond between aromatic and heterocycle
            if prod_mol.HasSubstructMatch(Chem.MolFromSmarts("[a][N]")):
                return True

        return False

    def dfs_traverse(node):
        nonlocal heterocycle_decoration_count

        if node["type"] == "reaction":
            rsmi = node["metadata"]["rsmi"]
            reactants = rsmi.split(">")[0].split(".")
            product = rsmi.split(">")[-1]

            # Check for heterocycle decoration
            if is_heterocycle_decoration(reactants, product):
                heterocycle_decoration_count += 1
                print(f"Heterocycle decoration detected: {reactants} -> {product}")

        # Recursively process children
        for child in node.get("children", []):
            dfs_traverse(child)

    # Start traversal from the root
    dfs_traverse(route)

    # Return True if the strategy is detected
    strategy_detected = heterocycle_decoration_count >= 2 and len(heterocycles_present) >= 2
    print(f"Heterocycle decoration strategy detected: {strategy_detected}")
    print(f"Heterocycles present: {heterocycles_present}")
    return strategy_detected
