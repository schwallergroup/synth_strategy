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
    Detects if the route contains a mid-stage N-arylation with a heterocycle.
    Mid-stage is defined as occurring in the middle of the synthesis.
    """
    found = False
    max_depth = 0

    def dfs_traverse(node, depth=0):
        nonlocal found, max_depth

        max_depth = max(max_depth, depth)

        if node["type"] == "reaction":
            if "rsmi" in node.get("metadata", {}):
                rsmi = node["metadata"]["rsmi"]
                reactants = rsmi.split(">")[0].split(".")
                product = rsmi.split(">")[-1]

                # Check for secondary amine in reactants
                amine_pattern = Chem.MolFromSmarts("[NH]")

                # Check for heterocycle patterns
                pyridine_pattern = Chem.MolFromSmarts("c1ccncc1")

                # Check for N-aryl bond in product
                n_aryl_pattern = Chem.MolFromSmarts("[N]-c:c")

                reactant_mols = [Chem.MolFromSmiles(r) for r in reactants if r]
                product_mol = Chem.MolFromSmiles(product) if product else None

                if product_mol and len(reactant_mols) >= 2:
                    has_amine = any(
                        mol.HasSubstructMatch(amine_pattern) for mol in reactant_mols if mol
                    )
                    has_heterocycle = any(
                        mol.HasSubstructMatch(pyridine_pattern) for mol in reactant_mols if mol
                    )

                    if (
                        has_amine
                        and has_heterocycle
                        and product_mol.HasSubstructMatch(n_aryl_pattern)
                    ):
                        # This is likely an N-arylation with a heterocycle
                        print(f"Found N-arylation with heterocycle at depth {depth}")

                        # Check if it's mid-stage (not first or last reaction)
                        if 0 < depth < max_depth:
                            found = True

        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    # First pass to determine max_depth
    def get_max_depth(node, depth=0):
        nonlocal max_depth
        max_depth = max(max_depth, depth)
        for child in node.get("children", []):
            get_max_depth(child, depth + 1)

    get_max_depth(route)
    mid_point = max_depth / 2

    # Second pass to find the pattern
    dfs_traverse(route)
    return found
