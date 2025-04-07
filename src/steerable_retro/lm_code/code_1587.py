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
    This function detects if the synthetic route involves modifications to
    multiple different ring systems (e.g., both benzene and pyrazole).
    """
    # Track modifications to different ring systems
    benzene_modified = False
    pyrazole_modified = False
    other_ring_modified = False

    def dfs_traverse(node, depth=0):
        nonlocal benzene_modified, pyrazole_modified, other_ring_modified

        if node["type"] == "reaction":
            if "rsmi" in node.get("metadata", {}):
                rsmi = node["metadata"]["rsmi"]
                reactants = rsmi.split(">")[0].split(".")
                product = rsmi.split(">")[-1]

                # Define patterns for different ring systems
                benzene_pattern = Chem.MolFromSmarts("c1ccccc1")
                pyrazole_pattern = Chem.MolFromSmarts("c1nn(C)cc1")  # Simplified N-methyl pyrazole

                product_mol = Chem.MolFromSmiles(product)
                if not product_mol:
                    return

                # Check for modifications by comparing with reactants
                for reactant in reactants:
                    reactant_mol = Chem.MolFromSmiles(reactant)
                    if not reactant_mol:
                        continue

                    # Check benzene modifications
                    if reactant_mol.HasSubstructMatch(
                        benzene_pattern
                    ) and product_mol.HasSubstructMatch(benzene_pattern):
                        # Compare benzene rings in reactant and product
                        if not Chem.MolToSmiles(reactant_mol) == Chem.MolToSmiles(product_mol):
                            benzene_modified = True
                            print(f"Benzene ring modified at depth {depth}")

                    # Check pyrazole modifications
                    if reactant_mol.HasSubstructMatch(
                        pyrazole_pattern
                    ) and product_mol.HasSubstructMatch(pyrazole_pattern):
                        # Compare pyrazole rings in reactant and product
                        if not Chem.MolToSmiles(reactant_mol) == Chem.MolToSmiles(product_mol):
                            pyrazole_modified = True
                            print(f"Pyrazole ring modified at depth {depth}")

                    # Check other ring modifications
                    # This is a simplified approach - in practice, you'd need more sophisticated
                    # ring detection and comparison
                    if "1" in reactant and "1" in product:  # Simple check for ring notation
                        if not (benzene_modified or pyrazole_modified):
                            other_ring_modified = True
                            print(f"Other ring system modified at depth {depth}")

        # Traverse children
        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    # Start traversal from root
    dfs_traverse(route)

    # Return True if at least two different ring systems were modified
    return (
        (benzene_modified and pyrazole_modified)
        or (benzene_modified and other_ring_modified)
        or (pyrazole_modified and other_ring_modified)
    )
