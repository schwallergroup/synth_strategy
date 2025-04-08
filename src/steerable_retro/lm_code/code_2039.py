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
    This function detects if the synthesis route involves coupling of a pyridine-containing
    fragment with another fragment.
    """
    result = False

    def dfs_traverse(node):
        nonlocal result

        if node["type"] == "reaction" and node.get("metadata", {}).get("rsmi"):
            rsmi = node["metadata"]["rsmi"]
            reactants = rsmi.split(">")[0].split(".")
            product = rsmi.split(">")[-1]

            # Check if there are at least 2 reactants (fragment coupling)
            if len(reactants) >= 2:
                # Create RDKit molecules
                try:
                    reactant_mols = [Chem.MolFromSmiles(r) for r in reactants]
                    prod_mol = Chem.MolFromSmiles(product)

                    # Check for pyridine in at least one reactant
                    pyridine_pattern = Chem.MolFromSmarts("c1ccncc1")
                    has_pyridine = any(
                        m and m.HasSubstructMatch(pyridine_pattern) for m in reactant_mols if m
                    )

                    # Check if the product is larger than each individual reactant
                    if has_pyridine and prod_mol:
                        prod_atoms = prod_mol.GetNumAtoms()
                        if all(m is None or m.GetNumAtoms() < prod_atoms for m in reactant_mols):
                            result = True
                            print("Pyridine fragment coupling detected")
                except:
                    pass

        # Traverse children
        for child in node.get("children", []):
            dfs_traverse(child)

    dfs_traverse(route)
    return result
