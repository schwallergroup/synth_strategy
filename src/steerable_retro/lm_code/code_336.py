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
    Detects a synthesis strategy that uses a bifunctional linker (like a dibromoalkane)
    to connect two different functional groups or heterocycles.
    """
    found_bifunctional_linker = False
    found_two_different_heterocycles = False

    # SMARTS patterns for common heterocycles
    heterocycle_patterns = {
        "morpholine": Chem.MolFromSmarts("[#7]1[#6][#6][#8][#6][#6]1"),
        "phthalimide": Chem.MolFromSmarts("c1ccccc1C(=O)NC(=O)"),
        "piperidine": Chem.MolFromSmarts("[#7]1[#6][#6][#6][#6][#6]1"),
        "pyrrolidine": Chem.MolFromSmarts("[#7]1[#6][#6][#6][#6]1"),
        "pyridine": Chem.MolFromSmarts("c1ccncc1"),
    }

    # Bifunctional linker patterns
    bifunctional_linker_patterns = [
        Chem.MolFromSmarts("[Br,Cl,I][#6]~[#6]~[#6]~[#6]~[#6][Br,Cl,I]"),  # Shorter chain
        Chem.MolFromSmarts("[Br,Cl,I][#6]~[#6]~[#6]~[#6]~[#6]~[#6]~[#6][Br,Cl,I]"),  # Longer chain
    ]

    heterocycles_found = set()

    def dfs_traverse(node):
        nonlocal found_bifunctional_linker, found_two_different_heterocycles, heterocycles_found

        if node["type"] == "reaction":
            if "metadata" in node and "rsmi" in node["metadata"]:
                rsmi = node["metadata"]["rsmi"]
                reactants = rsmi.split(">")[0].split(".")

                # Check for bifunctional linker in reactants
                for reactant in reactants:
                    try:
                        r_mol = Chem.MolFromSmiles(reactant)
                        if r_mol:
                            for pattern in bifunctional_linker_patterns:
                                if r_mol.HasSubstructMatch(pattern):
                                    found_bifunctional_linker = True
                                    print(f"Found bifunctional linker: {reactant}")
                    except:
                        continue

        elif node["type"] == "mol":
            # Check for heterocycles in molecules
            try:
                mol = Chem.MolFromSmiles(node["smiles"])
                if mol:
                    for hetero_name, pattern in heterocycle_patterns.items():
                        if mol.HasSubstructMatch(pattern):
                            heterocycles_found.add(hetero_name)
                            print(f"Found heterocycle: {hetero_name}")
            except:
                pass

        # Process children
        for child in node.get("children", []):
            dfs_traverse(child)

    # Start traversal
    dfs_traverse(route)

    # Check if we found at least two different heterocycles
    found_two_different_heterocycles = len(heterocycles_found) >= 2

    strategy_detected = found_bifunctional_linker and found_two_different_heterocycles
    print(f"Bifunctional linker strategy detected: {strategy_detected}")
    print(f"Heterocycles found: {heterocycles_found}")
    return strategy_detected
