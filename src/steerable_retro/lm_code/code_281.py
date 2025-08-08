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
    This function detects a synthetic strategy involving coupling of two heterocyclic fragments.
    """
    # Track if we found a coupling of heterocyclic fragments
    found_heterocycle_coupling = False

    def dfs_traverse(node):
        nonlocal found_heterocycle_coupling

        if node["type"] == "reaction":
            # Get reaction depth
            depth = int(node.get("metadata", {}).get("depth", -1))

            # Only consider final or near-final steps (depth 0 or 1)
            if depth <= 1:
                # Get reactants and product
                rsmi = node.get("metadata", {}).get("rsmi", "")
                if rsmi:
                    reactants_smiles = rsmi.split(">")[0].split(".")
                    product_smiles = rsmi.split(">")[-1]

                    # Define heterocycle patterns
                    thiazole_pattern = Chem.MolFromSmarts("c1scnc1")
                    indole_pattern = Chem.MolFromSmarts(
                        "c1cccc2c1[n]cc2"
                    )  # Matches both substituted and unsubstituted indoles

                    # Convert SMILES to molecules
                    reactant_mols = [Chem.MolFromSmiles(r) for r in reactants_smiles if r]

                    # Count heterocycles in reactants
                    heterocycle_count = 0
                    for mol in reactant_mols:
                        if mol:
                            if mol.HasSubstructMatch(thiazole_pattern):
                                heterocycle_count += 1
                                print(f"Found thiazole in reactant at depth {depth}")
                            if mol.HasSubstructMatch(indole_pattern):
                                heterocycle_count += 1
                                print(f"Found indole in reactant at depth {depth}")

                    # If we have at least 2 heterocycles in the reactants, it's a heterocycle coupling
                    if heterocycle_count >= 2:
                        found_heterocycle_coupling = True
                        print(f"Found heterocycle coupling at depth {depth}")

        # Traverse children
        for child in node.get("children", []):
            dfs_traverse(child)

    # Start traversal
    dfs_traverse(route)

    return found_heterocycle_coupling
